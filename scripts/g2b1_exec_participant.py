#!/usr/bin/env python3
"""External G2-B.1 executable-protocol participant (no deal.II linkage)."""
import argparse, os, socket, subprocess, sys, time
from pathlib import Path

MAGIC="G2B1_EXEC"; SCHEMA=1; LIMIT=64*1024*1024
HELLO=("dimension 3\nscalarBytes 8\nsignConvention R=Rs-fGamma;J=Js+JGamma\n"
       "forceUnits N\ntangentUnits N/m\ntimeIntegration BE,BDF2\nreplicatedMatrix 1\nend\n")

def fnv(data):
    value=1469598103934665603
    for byte in data: value=((value^byte)*1099511628211)&0xffffffffffffffff
    return value

def encode(kind,seq,payload,time_index=7,outer=3,operator=11,producer="simulatedParticipant"):
    payload=payload.encode(); pcheck=fnv(payload)
    prefix=f"{MAGIC} {SCHEMA} {kind} {seq} {len(payload)} {pcheck} {producer} dualConservative {time_index} {outer} {operator}"
    mcheck=fnv(prefix.encode()+b"\n"+payload)
    return f"{prefix} {mcheck}\n".encode()+payload,mcheck

def receive(sock):
    header=b""
    while not header.endswith(b"\n"):
        value=sock.recv(1)
        if not value: raise RuntimeError("truncated response header")
        header+=value
        if len(header)>4096: raise RuntimeError("response header too long")
    fields=header[:-1].decode().split()
    if len(fields)!=12 or fields[0]!=MAGIC or int(fields[1])!=SCHEMA: raise RuntimeError("bad response frame")
    length=int(fields[4]); payload=b""
    while len(payload)<length:
        value=sock.recv(length-len(payload))
        if not value: raise RuntimeError("truncated response payload")
        payload+=value
    if fnv(payload)!=int(fields[5]): raise RuntimeError("response payload checksum")
    prefix=" ".join(fields[:-1]).encode()
    if fnv(prefix+b"\n"+payload)!=int(fields[-1]): raise RuntimeError("response message checksum")
    return fields[2],int(fields[3]),payload.decode()

def send(sock,kind,seq,payload,**stamp):
    wire,checksum=encode(kind,seq,payload,**stamp); sock.sendall(wire); return checksum

def connect(path,timeout):
    deadline=time.monotonic()+timeout
    while True:
        try: client=socket.socket(socket.AF_UNIX,socket.SOCK_STREAM); client.connect(path); return client
        except (FileNotFoundError,ConnectionRefusedError):
            client.close()
            if time.monotonic()>=deadline: raise RuntimeError("connect timeout")
            time.sleep(.02)

def write_parameters(path,mesh,endpoint):
    path.write_text(f"""set mesh = {mesh}
set interface transfer = dualConservative
set protocol endpoint = {endpoint}
set protocol timeout seconds = 5
set rho = 970
set mu = 778200
set bulk modulus = 3631600
set delta t = 1e-7
set time integration = backwardEuler
set solid impedance = 0
set interface boundary = 4
set clamped boundary 1 = 1
set clamped boundary 2 = 2
set clamped boundary 3 = 3
set state input =
""")

def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--binary",required=True)
    parser.add_argument("--mesh",type=Path,required=True); parser.add_argument("--manifest",type=Path,required=True)
    parser.add_argument("--work",type=Path,required=True); parser.add_argument("--ranks",type=int,default=1)
    args=parser.parse_args(); args.work.mkdir(parents=True,exist_ok=True)
    endpoint=str(args.work/"protocol.sock"); prm=args.work/"protocol.prm"; write_parameters(prm,args.mesh,endpoint)
    command=(["mpirun","-np",str(args.ranks)] if args.ranks>1 else [])+[args.binary,str(prm)]
    log=(args.work/"solid.log").open("w"); process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
    try:
        client=connect(endpoint,8); client.settimeout(8)
        send(client,"Hello",1,HELLO)
        assert receive(client)[0]=="Capabilities"
        kind,_,dof_payload=receive(client); assert kind=="DofManifest"
        lines=dof_payload.splitlines(); dof_hash=lines[0].split()[1]
        entries=[line.split() for line in lines if line and line[0].isdigit()]
        free=next(values for values in entries if values[-1]=="0")
        free_entries=[values for values in entries if values[-1]=="0"]
        constrained=next(values for values in entries if values[-1]=="1")
        dof,component=int(free[0]),int(free[1])
        manifest=args.manifest.read_text(); graph=next(line.split()[1] for line in manifest.splitlines() if line.startswith("hashGraph "))
        operator_payload=(f"hashGraph {graph}\nhashWeights weights-v1\noperatorVersion 11\nzVersion 5\n"
          f"units H:1,W:m2\nsign fGamma=-HtWtf\nmanifestBytes {len(manifest.encode())}\n"+manifest)
        send(client,"OperatorManifest",2,operator_payload)
        kind,_,ready=receive(client); assert kind=="Ready" and dof_hash in ready
        common=f"zVersion 5\nhashGraph {graph}\nhashWeights weights-v1\ndofManifestHash {dof_hash}\n"
        force_payload=common+f"units N\nentries 1\n{dof} {component} 1\nend\n"
        tangent_payload=common+f"units N/m\nentries 1\n{dof} {dof} 100\nend\n"
        force_checksum=send(client,"ForceMessage",3,force_payload); assert receive(client)[0]=="Ack"
        tangent_checksum=send(client,"TangentMessage",4,tangent_payload); assert receive(client)[0]=="Ack"
        activation=(f"forceMessageChecksum {force_checksum}\ntangentMessageChecksum {tangent_checksum}\n"
                    f"hashGraph {graph}\nhashWeights weights-v1\nend\n")
        send(client,"ActivateCorrectorState",5,activation); assert receive(client)[0]=="Ack"
        kind,_,result=receive(client); assert kind=="AssemblyResult"
        values={line.split()[0]:line.split()[1] for line in result.splitlines() if len(line.split())==2}
        assert int(values["ranks"])==args.ranks
        assert float(values["rhsRepeatError"])==0 and float(values["jacobianRepeatError"])==0
        assert float(values["interfaceRhsNorm"])==1 and float(values["interfaceJacobianNorm"])==100
        assert int(values["interfaceJacobianEntries"])==1
        # Replay is ACK-idempotent and cannot reassemble.
        wire,_=encode("ActivateCorrectorState",5,activation); client.sendall(wire); assert receive(client)[0]=="Ack"
        send(client,"ClearProvisionalState",6,"end\n",outer=4); assert receive(client)[0]=="Ack"
        dof2,component2=int(free_entries[1][0]),int(free_entries[1][1])
        distributed_force=common+f"units N\nentries 2\n{dof} {component} 1\n{dof2} {component2} 2\nend\n"
        distributed_tangent=common+f"units N/m\nentries 2\n{dof} {dof} 100\n{dof2} {dof2} 200\nend\n"
        force_checksum=send(client,"ForceMessage",7,distributed_force,outer=4); assert receive(client)[0]=="Ack"
        tangent_checksum=send(client,"TangentMessage",8,distributed_tangent,outer=4); assert receive(client)[0]=="Ack"
        activation=(f"forceMessageChecksum {force_checksum}\ntangentMessageChecksum {tangent_checksum}\n"
                    f"hashGraph {graph}\nhashWeights weights-v1\nend\n")
        send(client,"ActivateCorrectorState",9,activation,outer=4); assert receive(client)[0]=="Ack"
        kind,_,distributed_result=receive(client); assert kind=="AssemblyResult"
        distributed_values={line.split()[0]:line.split()[1] for line in distributed_result.splitlines() if len(line.split())==2}
        assert abs(float(distributed_values["interfaceRhsNorm"])-5**.5)<1e-14
        assert abs(float(distributed_values["interfaceJacobianNorm"])-50000**.5)<1e-12

        send(client,"ClearProvisionalState",10,"end\n",outer=5); assert receive(client)[0]=="Ack"
        zero_force=common+f"units N\nentries 1\n{dof} {component} 0\nend\n"
        zero_tangent=common+f"units N/m\nentries 1\n{dof} {dof} 0\nend\n"
        fc=send(client,"ForceMessage",11,zero_force,outer=5); assert receive(client)[0]=="Ack"
        tc=send(client,"TangentMessage",12,zero_tangent,outer=5); assert receive(client)[0]=="Ack"
        activation=f"forceMessageChecksum {fc}\ntangentMessageChecksum {tc}\nhashGraph {graph}\nhashWeights weights-v1\nend\n"
        send(client,"ActivateCorrectorState",13,activation,outer=5); assert receive(client)[0]=="Ack"
        kind,_,zero_result=receive(client); assert kind=="AssemblyResult"
        zero_values={line.split()[0]:line.split()[1] for line in zero_result.splitlines() if len(line.split())==2}
        assert float(zero_values["interfaceRhsNorm"])==0 and float(zero_values["interfaceJacobianNorm"])==0

        send(client,"ClearProvisionalState",14,"end\n",outer=6); assert receive(client)[0]=="Ack"
        constrained_dof,constrained_component=int(constrained[0]),int(constrained[1])
        constrained_force=common+f"units N\nentries 1\n{constrained_dof} {constrained_component} 1\nend\n"
        constrained_tangent=common+f"units N/m\nentries 1\n{constrained_dof} {constrained_dof} 100\nend\n"
        fc=send(client,"ForceMessage",15,constrained_force,outer=6); assert receive(client)[0]=="Ack"
        tc=send(client,"TangentMessage",16,constrained_tangent,outer=6); assert receive(client)[0]=="Ack"
        activation=f"forceMessageChecksum {fc}\ntangentMessageChecksum {tc}\nhashGraph {graph}\nhashWeights weights-v1\nend\n"
        send(client,"ActivateCorrectorState",17,activation,outer=6); assert receive(client)[0]=="Ack"
        kind,_,constrained_result=receive(client); assert kind=="AssemblyResult"
        constrained_values={line.split()[0]:line.split()[1] for line in constrained_result.splitlines() if len(line.split())==2}
        assert float(constrained_values["interfaceRhsNorm"])==0 and float(constrained_values["interfaceJacobianNorm"])==0
        send(client,"Shutdown",18,"end\n",outer=6); assert receive(client)[0]=="Ack"; client.close()
        code=process.wait(timeout=10); assert code==0
        print(f"EXEC_PROTOCOL PASS ranks={args.ranks} dof={dof} rhsNorm={values['rhsNorm']} "
              f"repeat={values['rhsRepeatError']} Jrepeat={values['jacobianRepeatError']}")
    finally:
        if process.poll() is None: process.kill(); process.wait()
        log.close()
    return 0
if __name__=="__main__": raise SystemExit(main())
