#!/usr/bin/env python3
"""Validate the complete FILLETED_MASTER discrete interface without projection."""
import argparse, importlib.util, json, math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("geometry_audit",ROOT/"scripts/audit_g2a_interface_geometry.py")
audit=importlib.util.module_from_spec(spec); spec.loader.exec_module(audit)

def closest(point,quads,lower,upper,tol):
    mask=np.all((point>=lower-tol)&(point<=upper+tol),axis=1); indices=np.nonzero(mask)[0]
    if not len(indices): indices=np.arange(len(quads))
    best=None
    for index in indices:
        quad=quads[index]; _,area,normal,_=audit.polygon_stats(quad)
        signed=float((point-quad[0])@normal); projected=point-signed*normal
        # Convex planar polygon membership via oriented edge half-planes. This
        # avoids the cancellation suffered by a closest-triangle calculation
        # on 1.64-mm by sub-micrometre anisotropic master faces.
        side=[float(np.cross(quad[(i+1)%4]-quad[i],projected-quad[i])@normal) for i in range(4)]
        edge_scale=max(np.linalg.norm(quad[(i+1)%4]-quad[i]) for i in range(4))
        if not (min(side)>=-tol*edge_scale or max(side)<=tol*edge_scale): continue
        key=(abs(signed),int(index))
        if best is None or key<best[0]:
            e=point-projected; best=(key,{"distance":float(np.linalg.norm(e)),"q":projected,"e":e,"normal":normal,"signed":signed,"normal_distance":abs(signed),"tangential":0.0,"quad":int(index)})
    if best is not None: return best[1]
    local=audit.closest_on_quads(point,quads[indices]); local["quad"]=int(indices[local["quad"]]); return local

def main():
    parser=argparse.ArgumentParser(); parser.add_argument("bundle",type=Path); parser.add_argument("--contract",type=Path,default=ROOT/"geometry/interface-geometry-v1.json"); parser.add_argument("--json",type=Path,required=True); args=parser.parse_args()
    contract=json.loads(args.contract.read_text()); tol=contract["geometricTolerances"]["canonicalAbsolute"]
    _,nodes,surfaces,_=audit.gmsh(args.bundle/"solid.msh"); bid=contract["names"]["structuralInterfaceBoundaryId"]
    solid=np.array([[nodes[i] for i in conn] for tag,typ,conn in surfaces if tag==bid and typ==3]); slo=solid.min(1); shi=solid.max(1)
    poly=args.bundle/"convert/constant/polyMesh"; fp=audit.foam_points(poly/"points"); _,_,faces=audit.foam_patch(poly,contract["names"]["fluidInterfacePatch"]); fluid_faces=[fp[f] for f in faces]
    # closest_on_quads also handles a triangle represented with a repeated
    # fourth vertex; geometric area/normal calculations retain the true face.
    fluid=np.array([np.vstack((face,face[-1])) if len(face)==3 else face for face in fluid_faces]); flo=fluid.min(1); fhi=fluid.max(1)
    unique_ids=sorted(set(sum(faces,[]))); fs=[closest(fp[i],solid,slo,shi,tol) for i in unique_ids]
    solid_vertices=np.unique(solid.reshape(-1,3),axis=0); sf=[closest(p,fluid,flo,fhi,tol) for p in solid_vertices]
    fluid_area=0.; solid_area=0.; normal_error=0.; crossing=0; max_internal=0.; duplicate_faces=len(faces)-len({tuple(sorted(f)) for f in faces})
    for quad in solid: solid_area+=audit.polygon_stats(quad)[1]
    for face in fluid_faces:
        centre,area,nf,_=audit.polygon_stats(face); fluid_area+=area
        vertex_samples=list(face); interior_samples=[centre]
        if len(face)==4: interior_samples += [(face[0]+face[1]+face[2])/3,(face[0]+face[2]+face[3])/3]
        for point in vertex_samples:
            # Vertex membership is covered globally by fluidToSolid above;
            # shared master edges legitimately have more than one normal.
            closest(point,solid,slo,shi,tol)
        supports=[]
        for point in interior_samples:
            value=closest(point,solid,slo,shi,tol); max_internal=max(max_internal,value["distance"]); supports.append(value["quad"])
            ns=audit.polygon_stats(solid[value["quad"]])[2]; normal_error=max(normal_error,abs(float(nf@ns)+1.0))
        normals=[audit.polygon_stats(solid[i])[2] for i in sorted(set(supports))]
        if any(abs(float(normals[0]@normal))<1-1e-12 for normal in normals[1:]): crossing+=1
    fs_dist=[v["distance"] for v in fs]; sf_dist=[v["distance"] for v in sf]
    result={"classification":"PASS_FILLETED_MASTER","geometryContractHash":contract["geometryContractHash"],
      "fluidToSolid":audit.distribution(fs_dist),"solidToFluid":audit.distribution(sf_dist),
      "normalDistance":audit.distribution([v["normal_distance"] for v in fs]),"tangentialDistance":audit.distribution([v["tangential"] for v in fs]),
      "incompatiblePointCount":sum(v>tol for v in fs_dist),"fluidInterface":{"faces":len(faces),"points":len(unique_ids),"area":fluid_area},
      "solidMaster":{"faces":len(solid),"vertices":len(solid_vertices),"area":solid_area,"maximumQuadNonplanarity":max(abs(float((q[3]-q[0])@(np.cross(q[1]-q[0],q[2]-q[0])/np.linalg.norm(np.cross(q[1]-q[0],q[2]-q[0]))))) for q in solid)},
      "relativeAreaError":abs(fluid_area-solid_area)/solid_area,"maximumOppositeNormalError":normal_error,
      "maximumInteriorSampleDistance":max_internal,"nonCoplanarCrossingFaces":crossing,"duplicateFluidFaces":duplicate_faces,
      "gaps":sum(v>tol for v in fs_dist),"overlaps":duplicate_faces}
    required=(result["incompatiblePointCount"]==0 and max(sf_dist)<=tol and result["relativeAreaError"]<=contract["geometricTolerances"]["areaRelativeSerial"] and normal_error<=1e-12 and max_internal<=tol and crossing==0 and duplicate_faces==0)
    if not required: result["classification"]="BLOCKED_GEOMETRY"
    args.json.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n"); print(json.dumps(result,sort_keys=True)); return 0 if required else 1

if __name__=="__main__": raise SystemExit(main())
