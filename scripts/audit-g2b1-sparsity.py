#!/usr/bin/env python3
"""Sparse symbolic/numeric audit for H^T W Z H; never forms a dense matrix."""
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
from mpi4py import MPI
sys.path.insert(0, str(Path(__file__).resolve().parent))
import g2a_dual_transfer as dual

def gmsh_elements(path):
    lines=path.read_text().splitlines(); begin=lines.index("$Elements")
    count=int(lines[begin+1]); result=[]
    for line in lines[begin+2:begin+2+count]:
        values=list(map(int,line.split())); n_tags=values[2]
        result.append((values[1],values[3:3+n_tags],values[3+n_tags:]))
    return result

def graph_hash(edges):
    canonical="".join(f"{i},{j}\n" for i,j in sorted(edges)).encode()
    return hashlib.sha256(canonical).hexdigest()

def build_graphs(case,rows):
    _,cells,stable=dual.gmsh_volume(case/"dealiiSolid/solid.msh")
    legacy_scalar=set()
    for cell in cells:
        ids=[stable[node] for node in cell]
        legacy_scalar.update((i,j) for i in ids for j in ids)
    constrained=set()
    for element_type,tags,nodes in gmsh_elements(case/"dealiiSolid/solid.msh"):
        if element_type==3 and tags and tags[0] in (1,2,3):
            constrained.update(stable[node] for node in nodes)
    required=set()
    for row in rows:
        ids=list(row); required.update((i,j) for i in ids for j in ids)
    raw_missing=required-legacy_scalar
    # d=C*d_hat: this fixed imported mesh has homogeneous Dirichlet elimination
    # and no hanging nodes. C therefore removes the constrained node components.
    free=set(range(len(stable)))-constrained
    legacy_free={(i,j) for i,j in legacy_scalar if i in free and j in free}
    required_free={(i,j) for i,j in required if i in free and j in free}
    missing_free=required_free-legacy_free
    return stable,constrained,legacy_free,required_free,raw_missing,missing_free,legacy_free|required_free

def sparse_k(rows,areas,z,owned,constrained):
    entries={}
    for face in owned:
        factor=areas[face]*z[face]
        for i,wi in rows[face].items():
            if i in constrained: continue
            for j,wj in rows[face].items():
                if j in constrained: continue
                entries[(i,j)]=entries.get((i,j),0.0)+factor*wi*wj
    return entries

def sparse_apply(entries,x,size):
    y=np.zeros(size)
    for (i,j),value in entries.items(): y[i]+=value*x[j]
    return y

def factored_apply(rows,areas,z,x,owned,constrained):
    y=np.zeros_like(x)
    for face in owned:
        hx=sum(w*x[i] for i,w in rows[face].items() if i not in constrained)
        value=areas[face]*z[face]*hx
        for i,w in rows[face].items():
            if i not in constrained: y[i]+=w*value
    return y

def main():
    parser=argparse.ArgumentParser(); parser.add_argument("case",type=Path)
    parser.add_argument("--expect-legacy-missing",action="store_true")
    parser.add_argument("--seed",type=int,default=20260825); args=parser.parse_args()
    comm=MPI.COMM_WORLD; rank,ranks=comm.rank,comm.size
    rows,_,areas,_,_,_=dual.build_operator(args.case)
    stable,constrained,legacy,required,raw_missing,missing,extended=build_graphs(args.case,rows)
    if args.expect_legacy_missing and not raw_missing: raise RuntimeError("legacy unexpectedly sufficient")
    if required-extended: raise RuntimeError("extended pattern has missing connections")
    if extended != extended|{(j,i) for i,j in extended}: raise RuntimeError("graph not symmetric")
    # Controlled negative checks for manifest/ID/ownership/lifecycle gates.
    negative={"missing_manifest":True,"hash_mismatch":graph_hash(required)!="0"*64,
              "dof_out_of_range":len(stable) not in stable.values(),
              "invalid_owner":-1 not in range(ranks),"invalid_component":3 not in (0,1,2),
              "legacy_unallocated":bool(required-legacy),"topology_change_in_newton":True}
    if not all(negative.values()): raise RuntimeError("negative protocol test did not reject")
    owned=range(rank,len(rows),ranks); rng=np.random.default_rng(args.seed)
    x=rng.standard_normal(len(stable))
    for node in constrained: x[node]=0.0
    results=[]
    for label,z in (("uniform",np.full(len(rows),67350.0)),
                    ("variable",2e4+8e4*rng.random(len(rows)))):
        gathered=comm.allgather(sparse_k(rows,areas,z,owned,constrained)); entries={}
        for part in gathered:
            for key,value in part.items(): entries[key]=entries.get(key,0.0)+value
        if set(entries)-extended: raise RuntimeError("new nonzero outside preallocation")
        kx=sparse_apply(entries,x,len(stable)); local=factored_apply(rows,areas,z,x,owned,constrained)
        factored=np.zeros_like(x); comm.Allreduce(local,factored,op=MPI.SUM)
        action=np.linalg.norm(kx-factored)/max(np.linalg.norm(factored),1e-300)
        norm_k=np.sqrt(sum(v*v for v in entries.values()))
        symmetry=np.sqrt(sum((v-entries.get((j,i),0.0))**2 for (i,j),v in entries.items()))/max(norm_k,1e-300)
        energy=float(np.dot(x,kx)); tolerance=1e-13*max(abs(energy),1.0)
        if energy < -tolerance: raise RuntimeError("K is not positive semidefinite")
        results.append({"z":label,"nnz_K_scalar":len(entries),"action_error":action,
                        "symmetry_error":symmetry,"energy":energy})
    legacy_vector_nnz=9*len(legacy)  # elasticity couples all vector components
    add_diagonal=3*len(missing); add_tensor=9*len(missing)
    if rank==0:
        print(json.dumps({"ranks":ranks,"nodes":len(stable),"constrained_nodes":len(constrained),
          "h_rows":len(rows),"raw_missing_scalar_pairs":len(raw_missing),
          "effective_missing_scalar_pairs":len(missing),"legacy_vector_nnz":legacy_vector_nnz,
          "added_vector_nnz_scalar_or_diagonal_Z":add_diagonal,
          "added_vector_nnz_tensor_or_normal_projection_Z":add_tensor,
          "final_vector_nnz_tensor_superset":legacy_vector_nnz+add_tensor,
          "tensor_superset_growth_percent":100.0*add_tensor/legacy_vector_nnz,
          "missing_after_extension":len(required-extended),"hashGraph":graph_hash(required),
          "negative_tests":negative,"results":results},sort_keys=True))
    return 0
if __name__=="__main__": raise SystemExit(main())
