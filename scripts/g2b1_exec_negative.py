#!/usr/bin/env python3
"""Drive one controlled negative case through the real executable channel."""
import argparse, socket, subprocess, time
from pathlib import Path
import g2b1_exec_participant as proto

def expect_failure(process,client=None):
    if client:
        try:
            kind,_,_=proto.receive(client)
            if kind!="Nack": raise RuntimeError(f"expected Nack, got {kind}")
        except (BrokenPipeError,ConnectionError,TimeoutError,socket.timeout): pass
        client.close()
    code=process.wait(timeout=12)
    if code==0: raise RuntimeError("invalid protocol unexpectedly succeeded")

def launch(args,label):
    work=args.work/label; work.mkdir(parents=True,exist_ok=True)
    endpoint=str(work/"protocol.sock"); prm=work/"protocol.prm"
    proto.write_parameters(prm,args.mesh,endpoint)
    log=(work/"solid.log").open("w")
    process=subprocess.Popen([args.binary,str(prm)],stdout=log,stderr=subprocess.STDOUT)
    return process,log,endpoint

def handshake(client,manifest):
    proto.send(client,"Hello",1,proto.HELLO); assert proto.receive(client)[0]=="Capabilities"
    kind,_,dof_payload=proto.receive(client); assert kind=="DofManifest"
    lines=dof_payload.splitlines(); dof_hash=lines[0].split()[1]
    entries=[line.split() for line in lines if line and line[0].isdigit()]
    free=[values for values in entries if values[-1]=="0"]
    graph=next(line.split()[1] for line in manifest.splitlines() if line.startswith("hashGraph "))
    operator=(f"hashGraph {graph}\nhashWeights weights-v1\noperatorVersion 11\nzVersion 5\n"
      f"units H:1,W:m2\nsign fGamma=-HtWtf\nmanifestBytes {len(manifest.encode())}\n"+manifest)
    proto.send(client,"OperatorManifest",2,operator); assert proto.receive(client)[0]=="Ready"
    return dof_hash,graph,free

def main():
    parser=argparse.ArgumentParser(); parser.add_argument("case")
    parser.add_argument("--binary",required=True); parser.add_argument("--mesh",type=Path,required=True)
    parser.add_argument("--manifest",type=Path,required=True); parser.add_argument("--work",type=Path,required=True)
    args=parser.parse_args(); process,log,endpoint=launch(args,args.case)
    try:
        if args.case=="timeout": expect_failure(process); print("NEGATIVE PASS timeout"); return 0
        client=proto.connect(endpoint,8); client.settimeout(8)
        if args.case in {"bad_magic","bad_schema","unknown_type","oversize","bad_checksum","truncated"}:
            wire,_=proto.encode("Hello",1,proto.HELLO)
            if args.case=="bad_magic": wire=b"X"+wire[1:]
            elif args.case=="bad_schema": wire=wire.replace(b"G2B1_EXEC 1 ",b"G2B1_EXEC 2 ",1)
            elif args.case=="unknown_type": wire=wire.replace(b" Hello ",b" Unknown ",1)
            elif args.case=="oversize":
                fields=wire.split(b"\n",1)[0].split(); fields[4]=str(proto.LIMIT+1).encode()
                wire=b" ".join(fields)+b"\n"
            elif args.case=="bad_checksum": wire=wire[:-1]+bytes([wire[-1]^1])
            elif args.case=="truncated": wire=wire[:-4]
            client.sendall(wire)
            if args.case=="truncated": client.close(); client=None
            expect_failure(process,client); print("NEGATIVE PASS",args.case); return 0
        manifest=args.manifest.read_text(); dof_hash,graph,free=handshake(client,manifest)
        dof,component=int(free[0][0]),int(free[0][1]); other=int(free[-1][0])
        common=f"zVersion 5\nhashGraph {graph}\nhashWeights weights-v1\ndofManifestHash {dof_hash}\n"
        force=common+f"units N\nentries 1\n{dof} {component} 1\nend\n"
        tangent=common+f"units N/m\nentries 1\n{dof} {dof} 100\nend\n"
        if args.case=="nan": force=force.replace(" 1\nend"," nan\nend")
        if args.case=="inf": force=force.replace(" 1\nend"," inf\nend")
        if args.case=="duplicate_id": force=common+f"units N\nentries 2\n{dof} {component} 1\n{dof} {component} 2\nend\n"
        if args.case=="invalid_id": force=common+f"units N\nentries 1\n999999999 {component} 1\nend\n"
        if args.case=="invalid_tangent_id": tangent=common+"units N/m\nentries 1\n999999999 999999999 1\nend\n"
        if args.case=="force_pa": force=force.replace("units N\n","units Pa\n")
        if args.case=="tangent_units": tangent=tangent.replace("units N/m","units Pa")
        if args.case=="wrong_graph": force=force.replace(graph,"0"*64)
        if args.case=="wrong_weights": force=force.replace("weights-v1","wrong")
        if args.case=="wrong_dof_hash": force=force.replace(dof_hash,"wrong")
        if args.case=="outside_pattern": tangent=common+f"units N/m\nentries 2\n{dof} {other} 1\n{other} {dof} 1\nend\n"
        force_stamp={}; tangent_stamp={}
        if args.case=="stale": force_stamp["outer"]=2; tangent_stamp["outer"]=2
        if args.case=="future": force_stamp["outer"]=4; tangent_stamp["outer"]=4
        if args.case=="version_mismatch": tangent_stamp["operator"]=12
        if args.case=="tangent_only":
            tc=proto.send(client,"TangentMessage",3,tangent); assert proto.receive(client)[0]=="Ack"
            activation=(f"forceMessageChecksum 0\ntangentMessageChecksum {tc}\n"
                        f"hashGraph {graph}\nhashWeights weights-v1\nend\n")
            proto.send(client,"ActivateCorrectorState",4,activation)
            expect_failure(process,client); print("NEGATIVE PASS tangent_only"); return 0
        fc=proto.send(client,"ForceMessage",3,force,**force_stamp)
        kind,_,_=proto.receive(client)
        if kind=="Nack": expect_failure(process); print("NEGATIVE PASS",args.case); return 0
        if args.case=="sequence_collision":
            changed=force.replace(" 1\nend"," 2\nend"); wire,_=proto.encode("ForceMessage",3,changed)
            client.sendall(wire); expect_failure(process,client); print("NEGATIVE PASS sequence_collision"); return 0
        if args.case=="partial_disconnect": client.close(); expect_failure(process); print("NEGATIVE PASS partial_disconnect"); return 0
        if args.case!="force_only":
            tc=proto.send(client,"TangentMessage",4,tangent,**tangent_stamp)
            kind,_,_=proto.receive(client)
            if kind=="Nack": expect_failure(process); print("NEGATIVE PASS",args.case); return 0
        else: tc=0
        activation=(f"forceMessageChecksum {fc}\ntangentMessageChecksum {tc}\n"
                    f"hashGraph {graph}\nhashWeights weights-v1\nend\n")
        if args.case=="bad_activation": activation=activation.replace(str(fc),str(fc+1),1)
        proto.send(client,"ActivateCorrectorState",5,activation)
        expect_failure(process,client); print("NEGATIVE PASS",args.case)
    finally:
        if process.poll() is None: process.kill(); process.wait()
        log.close()
    return 0
if __name__=="__main__": raise SystemExit(main())
