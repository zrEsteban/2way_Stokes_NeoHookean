#!/usr/bin/env python3
"""Generate FILLETED_MASTER fluid/solid meshes from one versioned contract.

The deal.II boundary-4 Q1 faces are the discrete master.  The fluid boundary
is built from their exact end-section edges; no projection is used at runtime.
"""
import argparse, csv, hashlib, json, math, re, shutil, subprocess, tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT=ROOT/"geometry/interface-geometry-v1.json"
SOLID_TEMPLATE=ROOT/"cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS/dealiiSolid/roundedSolid.geo"
DEFAULT_CASE=ROOT/"cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS"

def canonical_hash(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),ensure_ascii=True).encode()).hexdigest()

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(1<<20),b""): h.update(block)
    return h.hexdigest()

def run(command,log):
    with log.open("w") as output:
        subprocess.run(command,check=True,stdout=output,stderr=subprocess.STDOUT)

def read_contract(path):
    data=json.loads(path.read_text()); claimed=data.pop("geometryContractHash")
    actual=canonical_hash(data)
    if claimed!=actual: raise RuntimeError(f"geometryContractHash mismatch: expected {actual}, got {claimed}")
    data["geometryContractHash"]=claimed
    if data["masterRepresentation"]!="FILLETED_MASTER" or data["projectionMode"]!="none":
        raise RuntimeError("this generator only implements FILLETED_MASTER without projection")
    return data

def parse_msh(path):
    lines=path.read_text().splitlines(); names={}
    i=lines.index("$PhysicalNames"); count=int(lines[i+1])
    for line in lines[i+2:i+2+count]:
        dim,tag,name=line.split(maxsplit=2); names[(int(dim),int(tag))]=name.strip('"')
    i=lines.index("$Nodes"); count=int(lines[i+1])
    nodes={int(v[0]):tuple(map(float,v[1:4])) for v in (line.split() for line in lines[i+2:i+2+count])}
    i=lines.index("$Elements"); count=int(lines[i+1]); surfaces=[]; volumes=[]
    for line in lines[i+2:i+2+count]:
        values=list(map(int,line.split())); eid,typ,nt=values[:3]; tags=values[3:3+nt]; conn=values[3+nt:]
        if typ==3 and tags: surfaces.append((eid,tags[0],conn))
        if typ==5: volumes.append((eid,conn))
    return names,nodes,surfaces,volumes

def master_cross_section(nodes,surfaces,boundary_id,length,tol=1e-14):
    # Boundary edges occurring once on the x=0 end of the extruded interface
    # form the open left-wall/roof/right-wall master polyline.
    edge_count={}
    for _,tag,conn in surfaces:
        if tag!=boundary_id: continue
        end=[nid for nid in conn if abs(nodes[nid][0])<=tol]
        if len(end)==2:
            edge=tuple(sorted(end)); edge_count[edge]=edge_count.get(edge,0)+1
    edges=[edge for edge,count in edge_count.items() if count==1]
    graph={}
    for a,b in edges: graph.setdefault(a,[]).append(b); graph.setdefault(b,[]).append(a)
    ends=[node for node,adj in graph.items() if len(adj)==1]
    if len(ends)!=2 or any(len(adj)>2 for adj in graph.values()): raise RuntimeError("boundary-4 end section is not one open polyline")
    start=min(ends,key=lambda nid:(nodes[nid][1],nodes[nid][2]))
    order=[start]; previous=None
    while True:
        nxt=[v for v in graph[order[-1]] if v!=previous]
        if not nxt: break
        previous,chosen=order[-1],nxt[0]; order.append(chosen)
    if len(order)!=len(graph): raise RuntimeError("incomplete master polyline")
    # The desired fluid loop runs glass left->right, then interface right->left.
    if nodes[order[0]][1] < nodes[order[-1]][1]: order=list(reversed(order))
    # A fluid face may cover several master faces only when they are strictly
    # coplanar. Remove intermediate collinear vertices but retain every Q1
    # normal transition (the chord junctions of both fillets).
    simplified=[order[0]]
    for i in range(1,len(order)-1):
        a=nodes[simplified[-1]]; b=nodes[order[i]]; c=nodes[order[i+1]]
        ab=[b[d]-a[d] for d in range(3)]; bc=[c[d]-b[d] for d in range(3)]
        cross=(ab[1]*bc[2]-ab[2]*bc[1],ab[2]*bc[0]-ab[0]*bc[2],ab[0]*bc[1]-ab[1]*bc[0])
        scale=math.sqrt(sum(v*v for v in ab))*math.sqrt(sum(v*v for v in bc))
        if math.sqrt(sum(v*v for v in cross))>1e-12*scale: simplified.append(order[i])
    simplified.append(order[-1])
    return simplified

def solid_wrapper(contract,output):
    d=contract["dimensions"]; f=contract["upperFillets"]; q=contract["discretization"]
    assignments={"L":d["length"],"w":d["width"],"h":d["height"],"r":f["radius"],
                 "t":d["sideWallThickness"],"e":d["roofThickness"],"nx":q["solidAxialCells"],
                 "lcNear":q["solidNearSize"],"lcFar":q["solidFarSize"],"allHex":1}
    output.write_text("// generated from interface-geometry-v1.json\n"+"\n".join(f"{k}={v:.17g};" for k,v in assignments.items())+f'\nInclude "{SOLID_TEMPLATE}";\n')

def fluid_geo(contract,nodes,order,output):
    d=contract["dimensions"]; q=contract["discretization"]; names=contract["names"]
    # Point 1/2 are glass endpoints. Remaining interface points are ordered
    # from right-bottom to left-bottom and reproduce every Q1 master segment.
    right=nodes[order[0]]; left=nodes[order[-1]]
    if abs(right[1]-d["width"])>1e-14 or abs(left[1])>1e-14: raise RuntimeError("unexpected interface endpoints")
    lines=["// generated from Q1 master surface; do not edit",'SetFactory("Built-in");',"Mesh.MshFileVersion=2.2;",f"Mesh.Algorithm={q['fluidMeshingAlgorithm']};","Mesh.RecombineAll=1;",f"Mesh.RecombinationAlgorithm={q['fluidRecombinationAlgorithm']};",f"Mesh.Smoothing={q['fluidSmoothingSteps']};"]
    radius=contract["upperFillets"]["radius"]
    point_coords=[left,(0.0,radius,0.0),(0.0,d["width"]-radius,0.0),right]+[nodes[nid] for nid in order[1:-1]]
    # loop: three glass segments left->right, then master chain right->left
    for index,p in enumerate(point_coords,1): lines.append(f"Point({index}) = {{{p[0]:.17g},{p[1]:.17g},{p[2]:.17g},{q['fluidInteriorSize']:.17g}}};")
    polygon=list(range(1,len(point_coords)+1))+[1]
    for line_id,(a,b) in enumerate(zip(polygon,polygon[1:]),1): lines.append(f"Line({line_id})={{{a},{b}}};")
    curve_ids=list(range(1,len(polygon))); lines.append("Curve Loop(1)={"+",".join(map(str,curve_ids))+"};")
    if len(point_coords)!=10: raise RuntimeError("FILLETED_MASTER expects two Q1 chords per upper fillet")
    width=q["fluidWidthCells"]; height=q["fluidHeightCells"]
    if width<4 or height<3 or q["filletChordCells"]!=1: raise RuntimeError("invalid canonical transfinite discretization")
    # Four logical sides: glass; right wall+first chord; second chord+roof+
    # second chord; first chord+left wall. Composite sides have matching totals.
    progression=q["fluidSideProgression"]
    lines.extend(["Transfinite Curve{1,3,5,6,8,9}=2;",f"Transfinite Curve{{2,7}}={width-1};",f"Transfinite Curve{{4}}={height} Using Progression {progression:.17g};",f"Transfinite Curve{{10}}={height} Using Progression {1/progression:.17g};",
                  "Plane Surface(1)={1};","Transfinite Surface{1}={1,4,6,9};","Recombine Surface{1};",f"out[]=Extrude {{{d['length']:.17g},0,0}} {{ Surface{{1}}; Layers{{{q['fluidAxialCells']}}}; Recombine; }};",
                  'Physical Volume("fluid")={out[1]};','Physical Surface("inlet")={1};','Physical Surface("outlet")={out[0]};','Physical Surface("glass")={out[2],out[3],out[4]};',
                  'Physical Surface("interface")={'+",".join(f"out[{i}]" for i in range(5,len(curve_ids)+2))+'};'])
    output.write_text("\n".join(lines)+"\n")

def foam_points(path):
    text=path.read_text(); return [tuple(map(float,m)) for m in re.findall(r"\(([-+0-9.eE]+) ([-+0-9.eE]+) ([-+0-9.eE]+)\)",text)]

def foam_faces(path):
    return [list(map(int,v.split())) for v in re.findall(r"\d+\(([^()]*)\)",path.read_text())]

def foam_patch(mesh,name):
    text=(mesh/"boundary").read_text(); match=re.search(rf"\b{name}\s*\{{.*?nFaces\s+(\d+);.*?startFace\s+(\d+);",text,re.S)
    if not match: raise RuntimeError(f"missing OpenFOAM patch {name}")
    count,start=map(int,match.groups()); return start,foam_faces(mesh/"faces")[start:start+count]

def write_queries_and_metadata(poly_mesh,contract,nodes,surfaces,output):
    points=foam_points(poly_mesh/"points"); start,faces=foam_patch(poly_mesh,contract["names"]["fluidInterfacePatch"])
    ids=[]; seen=set()
    for face in faces:
        for point in face:
            if point not in seen: seen.add(point); ids.append(point)
    with (output/"robin-query.csv").open("w",newline="") as stream:
        writer=csv.writer(stream); writer.writerow(["x","y","z"]); writer.writerows(points[i] for i in ids)
    master=[]
    for eid,tag,conn in surfaces:
        if tag==contract["names"]["structuralInterfaceBoundaryId"]: master.append((eid,conn,[nodes[n] for n in conn]))
    # Metadata is diagnostic only. Exact membership is checked against every
    # planar Q1 master face and the deterministic lowest element ID is used.
    rows=[]; tol=contract["geometricTolerances"]["canonicalAbsolute"]
    for local,pid in enumerate(ids):
        p=points[pid]; candidates=[]
        for eid,conn,coords in master:
            lo=[min(v[d] for v in coords)-tol for d in range(3)]; hi=[max(v[d] for v in coords)+tol for d in range(3)]
            if all(lo[d]<=p[d]<=hi[d] for d in range(3)):
                candidates.append((eid,conn))
        if not candidates: raise RuntimeError(f"generated fluid point {pid} is not on Q1 master bbox")
        eid,conn=min(candidates)
        rows.append({"queryPointId":local,"globalFluidPointId":pid,"masterElementId":eid,"masterFaceIndex":4,
                     "physicalCoordinate":p,"geometryVersion":contract["geometryVersion"]})
    (output/"interface-generation-metadata.json").write_text(json.dumps(rows,indent=2,sort_keys=True)+"\n")
    return len(points),len(faces),len(ids)

def semantic_mesh_hash(poly_mesh):
    points=foam_points(poly_mesh/"points"); faces=foam_faces(poly_mesh/"faces")
    boundary=re.sub(r"//.*","",(poly_mesh/"boundary").read_text())
    return canonical_hash({"points":points,"faces":faces,"boundary":" ".join(boundary.split())})

def generate(contract_path,case,output,install=False):
    contract=read_contract(contract_path); output.mkdir(parents=True,exist_ok=False)
    logs=output/"logs"; logs.mkdir(); wrapper=output/"solid.generated.geo"; solid=output/"solid.msh"
    solid_wrapper(contract,wrapper); run(["gmsh",str(wrapper),"-3","-format","msh2","-o",str(solid),"-v","2"],logs/"gmsh-solid.log")
    names,nodes,surfaces,volumes=parse_msh(solid); bid=contract["names"]["structuralInterfaceBoundaryId"]
    if names.get((2,bid))!=contract["names"]["structuralInterfacePhysicalName"]: raise RuntimeError("structural boundary identity mismatch")
    order=master_cross_section(nodes,surfaces,bid,contract["dimensions"]["length"])
    fluid_source=output/"fluid.generated.geo"; fluid_msh=output/"fluid.msh"; fluid_geo(contract,nodes,order,fluid_source)
    run(["gmsh",str(fluid_source),"-3","-format","msh2","-o",str(fluid_msh),"-v","2"],logs/"gmsh-fluid.log")
    _,fluid_nodes,_,fluid_volumes=parse_msh(fluid_msh)
    convert=output/"convert"; (convert/"system").mkdir(parents=True)
    control=(case/"system/controlDict").read_text()
    control=re.sub(r"writePrecision\s+\d+\s*;","writePrecision 17;",control)
    (convert/"system/controlDict").write_text(control)
    for dictionary in ("fvSchemes","fvSolution"):
        source=case/"system/fluid"/dictionary
        if source.exists(): shutil.copy2(source,convert/"system"/dictionary)
    run(["gmshToFoam","-case",str(convert),str(fluid_msh)],logs/"gmshToFoam-fluid.log")
    poly=convert/"constant/polyMesh"; point_count,interface_faces,interface_points=write_queries_and_metadata(poly,contract,nodes,surfaces,output)
    manifest={"schemaVersion":1,"geometryContractHash":contract["geometryContractHash"],"geometryVersion":contract["geometryVersion"],
      "tools":{"gmsh":subprocess.run(["gmsh","--version"],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=True).stdout.strip()},
      "solid":{"cells":len(volumes),"points":len(nodes),"interfaceFaces":sum(tag==bid for _,tag,_ in surfaces),"sha256":sha(solid)},
      "fluid":{"cells":len(fluid_volumes),"points":point_count,"gmshPoints":len(fluid_nodes),"interfaceFaces":interface_faces,"interfacePoints":interface_points,"mshSha256":sha(fluid_msh),"polyMeshSemanticHash":semantic_mesh_hash(poly)},
      "inputs":{"contractSha256":sha(contract_path),"solidTemplateSha256":sha(SOLID_TEMPLATE),"generatorSha256":sha(Path(__file__))},"masterCrossSectionPoints":len(order)}
    manifest["manifestHash"]=canonical_hash(manifest); (output/"generation-manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    if install:
        fluid_dest=case/"constant/fluid/polyMesh"; shutil.rmtree(fluid_dest,ignore_errors=True); shutil.copytree(poly,fluid_dest)
        shutil.copy2(solid,case/"dealiiSolid/solid.msh"); shutil.copy2(output/"robin-query.csv",case/"dealiiSolid/robin-query.csv")
        shutil.copy2(output/"interface-generation-metadata.json",case/"dealiiSolid/interface-generation-metadata.json")
    print(json.dumps(manifest,sort_keys=True))

if __name__=="__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("--contract",type=Path,default=DEFAULT_CONTRACT); parser.add_argument("--case",type=Path,default=DEFAULT_CASE); parser.add_argument("--output",type=Path,required=True); parser.add_argument("--install",action="store_true")
    args=parser.parse_args(); generate(args.contract.resolve(),args.case.resolve(),args.output.resolve(),args.install)
