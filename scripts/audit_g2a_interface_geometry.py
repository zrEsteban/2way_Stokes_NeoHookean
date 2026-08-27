#!/usr/bin/env python3
"""Read-only G2-A geometry audit; never constructs or changes H_ps."""
import argparse, csv, hashlib, json, math, os, re
from pathlib import Path
import numpy as np

def foam_points(path):
    text=path.read_text()
    return np.array([tuple(map(float,m)) for m in re.findall(
        r"\(([-+0-9.eE]+) ([-+0-9.eE]+) ([-+0-9.eE]+)\)",text)])

def foam_faces(path):
    text=path.read_text()
    return [list(map(int,v.split())) for v in re.findall(r"\d+\(([^()]*)\)",text)]

def foam_patch(mesh,name):
    boundary=(mesh/"boundary").read_text()
    names=[value for value in re.findall(r"^\s*([A-Za-z0-9_]+)\s*\{",boundary,re.M) if value!="FoamFile"]
    match=re.search(rf"\b{re.escape(name)}\s*\{{.*?nFaces\s+(\d+);.*?startFace\s+(\d+);",boundary,re.S)
    if not match: raise RuntimeError(f"patch {name} not found")
    n,start=map(int,match.groups()); all_faces=foam_faces(mesh/"faces")
    return names.index(name),start,all_faces[start:start+n]

def gmsh(path):
    lines=path.read_text().splitlines(); names={}
    p=lines.index("$PhysicalNames"); n=int(lines[p+1])
    for line in lines[p+2:p+2+n]:
        dim,tag,name=line.split(maxsplit=2); names[(int(dim),int(tag))]=name.strip('"')
    p=lines.index("$Nodes"); n=int(lines[p+1]); nodes={int(a[0]):np.array(list(map(float,a[1:4]))) for a in (x.split() for x in lines[p+2:p+2+n])}
    p=lines.index("$Elements"); n=int(lines[p+1]); surfaces=[]; volume_types=set()
    for line in lines[p+2:p+2+n]:
        a=list(map(int,line.split())); typ,nt=a[1],a[2]; tags=a[3:3+nt]; conn=a[3+nt:]
        if typ in (2,3) and tags: surfaces.append((tags[0],typ,conn))
        if typ in (4,5,6,11,12,17): volume_types.add(typ)
    return names,nodes,surfaces,sorted(volume_types)

def triangle_closest(p,a,b,c):
    # Ericson, Real-Time Collision Detection.
    ab=b-a; ac=c-a; ap=p-a; d1=ab@ap; d2=ac@ap
    if d1<=0 and d2<=0:return a,(1.,0.,0.)
    bp=p-b; d3=ab@bp; d4=ac@bp
    if d3>=0 and d4<=d3:return b,(0.,1.,0.)
    vc=d1*d4-d3*d2
    if vc<=0 and d1>=0 and d3<=0:
        v=d1/(d1-d3); return a+v*ab,(1-v,v,0.)
    cp=p-c; d5=ab@cp; d6=ac@cp
    if d6>=0 and d5<=d6:return c,(0.,0.,1.)
    vb=d5*d2-d1*d6
    if vb<=0 and d2>=0 and d6<=0:
        w=d2/(d2-d6); return a+w*ac,(1-w,0.,w)
    va=d3*d6-d5*d4
    if va<=0 and (d4-d3)>=0 and (d5-d6)>=0:
        w=(d4-d3)/((d4-d3)+(d5-d6)); return b+w*(c-b),(0.,1-w,w)
    denom=1/(va+vb+vc); v=vb*denom; w=vc*denom
    return a+ab*v+ac*w,(1-v-w,v,w)

def polygon_stats(x):
    edges=np.roll(x,-1,axis=0)-x; lengths=np.linalg.norm(edges,axis=1)
    centre=x.mean(axis=0); av=np.zeros(3)
    for i in range(len(x)): av+=np.cross(x[i]-centre,x[(i+1)%len(x)]-centre)/2
    area=np.linalg.norm(av); normal=av/max(area,1e-300)
    return centre,area,normal,lengths

def point_segment_distance(point,a,b):
    direction=b-a; denominator=float(direction@direction)
    parameter=0. if denominator==0 else min(1.,max(0.,float((point-a)@direction)/denominator))
    return float(np.linalg.norm(point-(a+parameter*direction)))

def closest_on_quads(point,quads):
    best=None
    for qi,q in enumerate(quads):
        for tri,ids in enumerate(((0,1,2),(0,2,3))):
            closest,bary=triangle_closest(point,q[ids[0]],q[ids[1]],q[ids[2]])
            d=np.linalg.norm(point-closest)
            key=(d,qi,tri)
            if best is None or key<best[0]: best=(key,closest,bary,ids)
    _,q,bary,ids=best; qi=best[0][1]; tri=best[0][2]
    centre,area,normal,lengths=polygon_stats(quads[qi]); e=point-q
    signed=float(e@normal); tang=e-signed*normal
    location="interior" if min(bary)>1e-10 else ("vertex" if sum(v<1e-10 for v in bary)>=2 else "edge")
    # Diagnostic bilinear coordinates for the selected split.
    if tri==0: uv=(bary[1]+bary[2],bary[2])
    else: uv=(bary[1],bary[1]+bary[2])
    return dict(distance=float(np.linalg.norm(e)),q=q,e=e,normal=normal,signed=signed,
                normal_distance=abs(signed),tangential=float(np.linalg.norm(tang)),
                quad=qi,triangle=tri,uv=uv,location=location,area=area,lengths=lengths)

def distribution(values):
    a=np.asarray(values,float)
    return {"min":float(a.min()),"max":float(a.max()),"mean":float(a.mean()),
            "rms":float(np.sqrt(np.mean(a*a))),**{f"p{p}":float(np.percentile(a,p)) for p in (50,90,95,99)}}

def rigid_fit(source,target,similarity=False):
    cs=source.mean(0); ct=target.mean(0); xs=source-cs; xt=target-ct
    u,s,vt=np.linalg.svd(xs.T@xt); r=vt.T@u.T
    if np.linalg.det(r)<0: vt[-1]*=-1; r=vt.T@u.T
    scale=float(s.sum()/np.sum(xs*xs)) if similarity else 1.
    t=ct-scale*r@cs; fit=(scale*(r@source.T)).T+t; residual=np.linalg.norm(fit-target,axis=1)
    return {"scale":scale,"rotation":r.tolist(),"translation":t.tolist(),"residual":distribution(residual),"condition":float(s[0]/max(s[-1],1e-300))}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("case",type=Path); ap.add_argument("--json",type=Path,required=True); ap.add_argument("--csv",type=Path,required=True)
    args=ap.parse_args(); case=args.case; mesh=case/"constant/fluid/polyMesh"
    fp=foam_points(mesh/"points"); patch_id,start,faces=foam_patch(mesh,"interface")
    patch_point_ids=sorted(set(sum(faces,[]))); sorted_fluid=fp[patch_point_ids]
    patch_local_order=[]; seen_points=set()
    for face in faces:
        for point_id in face:
            if point_id not in seen_points: seen_points.add(point_id); patch_local_order.append(point_id)
    fixture=np.loadtxt(case/"dealiiSolid/robin-query.csv",delimiter=",",skiprows=1)
    coordinate_to_global={tuple(fp[i]):i for i in patch_point_ids}
    fixture_ids=[coordinate_to_global.get(tuple(p)) for p in fixture]
    fixture_set_match=(None not in fixture_ids and len(set(fixture_ids))==len(patch_point_ids)
                       and set(fixture_ids)==set(patch_point_ids))
    fixture_order_match=fixture_set_match and fixture_ids==patch_point_ids
    fixture_patch_local_order_match=fixture_set_match and fixture_ids==patch_local_order
    if not fixture_set_match: raise RuntimeError("query fixture is not an exact permutation of OpenFOAM patch points")
    point_ids=fixture_ids; fluid=fixture.copy()
    names,nodes,surfaces,volume_types=gmsh(case/"dealiiSolid/solid.msh")
    boundary={}
    for tag,typ,conn in surfaces:
        x=np.array([nodes[i] for i in conn]); boundary.setdefault(tag,[]).append(x)
    solid_quads=np.array([q for q in boundary[4] if len(q)==4])
    if len(solid_quads)!=len(boundary[4]): raise RuntimeError("boundary 4 is not all Q1 quadrilaterals")
    face_nonplanarity=[]
    for quad in solid_quads:
        normal=np.cross(quad[1]-quad[0],quad[2]-quad[0]); normal/=np.linalg.norm(normal)
        face_nonplanarity.append(abs(float((quad[3]-quad[0])@normal)))

    solid_lower=solid_quads.min(axis=1); solid_upper=solid_quads.max(axis=1)
    fs=[]; accelerated_equivalent=0; accelerated_same_face=0
    for i,p in enumerate(fluid):
        value=closest_on_quads(p,solid_quads); value["query_id"]=i; value["global_fluid_point_id"]=point_ids[i]
        box_delta=np.maximum(np.maximum(solid_lower-p,0),p-solid_upper)
        box_distance=np.linalg.norm(box_delta,axis=1); accelerated=int(np.argmin(box_distance))
        accelerated_value=closest_on_quads(p,solid_quads[accelerated:accelerated+1])
        value["acceleratedFace"]=accelerated; value["acceleratedDistance"]=accelerated_value["distance"]
        value["acceleratedEquivalent"]=abs(accelerated_value["distance"]-value["distance"])<=1e-14
        accelerated_equivalent+=value["acceleratedEquivalent"]; accelerated_same_face+=accelerated==value["quad"]
        fs.append(value)
    solid_vertices=np.unique(solid_quads.reshape(-1,3),axis=0)
    fluid_quads=np.array([fp[f] for f in faces])
    sf=[closest_on_quads(p,fluid_quads) for p in solid_vertices]

    boundary_table=[]
    for tag,polys in sorted(boundary.items()):
        areas=[]; centres=[]; normals=[]; vertices=[]; edges=[]
        for poly in polys:
            c,a,n,l=polygon_stats(poly); areas.append(a); centres.append(c); normals.append(n*a); vertices.extend(poly); edges.extend(l)
        v=np.array(vertices); total=sum(areas); centroid=np.average(np.array(centres),axis=0,weights=areas)
        normal=np.sum(normals,axis=0); normal/=max(np.linalg.norm(normal),1e-300)
        boundary_table.append({"boundaryId":tag,"physicalName":names.get((2,tag),"unknown"),"faceCount":len(polys),"area":total,
          "bboxMin":v.min(0).tolist(),"bboxMax":v.max(0).tolist(),"centroid":centroid.tolist(),"meanNormal":normal.tolist(),
          "minEdge":min(edges),"maxEdge":max(edges)})

    # Patch aggregate and provenance.
    fstats=[polygon_stats(q) for q in fluid_quads]; farea=sum(v[1] for v in fstats)
    fcent=np.average(np.array([v[0] for v in fstats]),axis=0,weights=[v[1] for v in fstats])
    fnormal=sum((v[2]*v[1] for v in fstats),np.zeros(3)); fnormal/=np.linalg.norm(fnormal)
    sqstats=[polygon_stats(q) for q in solid_quads]; sarea=sum(v[1] for v in sqstats)
    scent=np.average(np.array([v[0] for v in sqstats]),axis=0,weights=[v[1] for v in sqstats])

    nearest=np.array([v["q"] for v in fs]); before=np.array([v["distance"] for v in fs])
    translation=(nearest-fluid).mean(0); translated=np.linalg.norm(fluid+translation-nearest,axis=1)
    signed=np.array([v["signed"] for v in fs]); normal_offset=float(signed.mean())
    axis_scale=[]
    for d in range(3):
        A=np.column_stack((fluid[:,d],np.ones(len(fluid)))); axis_scale.append(np.linalg.lstsq(A,nearest[:,d],rcond=None)[0].tolist())
    axis_fit=np.column_stack([axis_scale[d][0]*fluid[:,d]+axis_scale[d][1] for d in range(3)])

    rng=np.random.default_rng(20260827); forces=rng.standard_normal((len(fluid),3)); omegas=np.eye(3)
    moment_before=np.cross(fluid,forces).sum(0); moment_after=np.cross(nearest,forces).sum(0)
    moment_error=np.linalg.norm(moment_before-moment_after)/max(np.linalg.norm(moment_before),1e-300)
    rotation_errors=[]; work_errors=[]
    for omega in omegas:
        df=np.cross(np.broadcast_to(omega,fluid.shape),fluid); dq=np.cross(np.broadcast_to(omega,nearest.shape),nearest)
        rotation_errors.append(float(np.max(np.linalg.norm(df-dq,axis=1))))
        wf=float(np.sum(df*forces)); wq=float(np.sum(dq*forces)); work_errors.append(abs(wf-wq)/max(abs(wf),abs(wq),1e-300))

    rank=int(os.environ.get("OMPI_COMM_WORLD_RANK",os.environ.get("PMI_RANK","0")))
    with (args.csv.open("w",newline="") if rank==0 else open(os.devnull,"w",newline="")) as out:
        writer=csv.writer(out); writer.writerow(["queryPointId","globalFluidPointId","x","y","z","qx","qy","qz","ex","ey","ez","distance","signedNormalDistance","normalDistance","tangentialDistance","solidFaceIndex","acceleratedFaceIndex","acceleratedDistance","acceleratedEquivalent","triangle","u","v","location","region","hDiameter","hMinEdge","hMaxEdge","hArea","aspect","dOverDiameter","dOverMinEdge","dOverArea","dOverInterfaceLength","dOverCoordinateEpsilon"])
        L=np.linalg.norm(fluid.max(0)-fluid.min(0))
        for v in fs:
            hdiam=max(np.linalg.norm(solid_quads[v["quad"]][i]-solid_quads[v["quad"]][j]) for i in range(4) for j in range(i)); hmin=min(v["lengths"]); hmax=max(v["lengths"]); harea=math.sqrt(v["area"])
            p=fluid[v["query_id"]]
            if abs(p[2]-3.3e-5)<1e-15 and abs(p[1])<1e-15: region="upper-left-sharp-edge"
            elif abs(p[2]-3.3e-5)<1e-15 and abs(p[1]-4.92e-4)<1e-15: region="upper-right-sharp-edge"
            elif abs(p[2]-3.3e-5)<1e-15: region="roof-interior"
            elif abs(p[1])<1e-15 or abs(p[1]-4.92e-4)<1e-15: region="side-interior"
            else: region="other"
            coordinate_epsilon=np.finfo(float).eps*max(np.max(np.abs(fluid)),1e-300)
            writer.writerow([v["query_id"],v["global_fluid_point_id"],*p,*v["q"],*v["e"],v["distance"],v["signed"],v["normal_distance"],v["tangential"],v["quad"],v["acceleratedFace"],v["acceleratedDistance"],int(v["acceleratedEquivalent"]),v["triangle"],*v["uv"],v["location"],region,hdiam,hmin,hmax,harea,hmax/hmin,v["distance"]/hdiam,v["distance"]/hmin,v["distance"]/harea,v["distance"]/L,v["distance"]/coordinate_epsilon])

    result={
      "classification":"CAD_MISMATCH",
      "provenance":{"fluidPatch":"interface","fluidPatchId":patch_id,"startFace":start,"faceCount":len(faces),"uniquePointCount":len(point_ids),"fixtureExactSetMatch":fixture_set_match,"fixtureGlobalOrderMatch":fixture_order_match,"fixturePatchLocalOrderMatch":fixture_patch_local_order_match,"coordinateFrame":"reference","units":"m","query20GlobalFluidPointId":point_ids[20],"query20Faces":[start+i for i,f in enumerate(faces) if point_ids[20] in f]},
      "fluidSurface":{"bboxMin":fluid.min(0).tolist(),"bboxMax":fluid.max(0).tolist(),"area":farea,"centroid":fcent.tolist(),"meanNormal":fnormal.tolist()},
      "solidBoundary4":{"bboxMin":solid_vertices.min(0).tolist(),"bboxMax":solid_vertices.max(0).tolist(),"area":sarea,"centroid":scent.tolist(),"faceCount":len(solid_quads),"maximumQuadNonplanarity":max(face_nonplanarity)},
      "boundaryTable":boundary_table,"gmshVolumeElementTypes":volume_types,"mapping":"Q1 over type-5 hexahedra; surface type-3 quads",
      "fluidToSolid":distribution(before),"solidToFluid":distribution([v["distance"] for v in sf]),
      "normalDistance":distribution([v["normal_distance"] for v in fs]),"tangentialDistance":distribution([v["tangential"] for v in fs]),"signedNormalDistance":distribution(signed),
      "searchComparison":{"points":len(fs),"sameFace":accelerated_same_face,"equivalentDistance":accelerated_equivalent,"differentDistance":len(fs)-accelerated_equivalent},
      "outsideTolerance":{"1e-12":sum(v["distance"]>1e-12 for v in fs),"1e-10":sum(v["distance"]>1e-10 for v in fs),"1e-8":sum(v["distance"]>1e-8 for v in fs)},
      "query20":{k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in fs[20].items() if k not in ("lengths",)},
      "diagnosticFits":{"translation":{"value":translation.tolist(),"residual":distribution(translated)},"rigid":rigid_fit(fluid,nearest),"similarity":rigid_fit(fluid,nearest,True),"axisScaleAndOffset":axis_scale,"axisResidual":distribution(np.linalg.norm(axis_fit-nearest,axis=1)),"constantNormalOffset":{"value":normal_offset,"residual":distribution(np.abs(signed-normal_offset))}},
      "projectionImpact":{"randomForceMomentRelativeError":moment_error,"unitAxisRotationMaxAbs":rotation_errors,"unitAxisRotationWorkRelativeError":work_errors},
      "nominalSources":{"fluid":"blockMesh sharp rectangular channel interface","solid":"OpenCASCADE roundedSolid.geo with r=7.5e-7 m circular upper fillets","commonSurface":False},
    }
    query20_quad=solid_quads[result["query20"]["quad"]]
    result["query20"]["distanceToFaceEdges"]=[point_segment_distance(result["query20"]["q"],query20_quad[i],query20_quad[(i+1)%4]) for i in range(4)]
    canonical=json.dumps(result,sort_keys=True,separators=(",",":")); result["auditHash"]=hashlib.sha256(canonical.encode()).hexdigest()
    if rank==0: args.json.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"rank":rank,"classification":result["classification"],"auditHash":result["auditHash"],"fluidToSolid":result["fluidToSolid"],"query20":result["query20"],"fixtureExactSetMatch":fixture_set_match,"fixtureGlobalOrderMatch":fixture_order_match,"fixturePatchLocalOrderMatch":fixture_patch_local_order_match},sort_keys=True))

if __name__=="__main__": raise SystemExit(main())
