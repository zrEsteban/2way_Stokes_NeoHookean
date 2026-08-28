#!/usr/bin/env python3
"""Summarize checkMesh fields and spatially map low-determinant cells."""
import argparse, csv, json, math, re
from pathlib import Path

THRESHOLDS=(1e-5,5e-5,1e-4,5e-4,1e-3,2e-3,5e-3,1e-2)

def values(path, vector=False):
    text=path.read_text(); match=re.search(r"internalField\s+nonuniform\s+List<[^>]+>\s+\d+\s*\((.*?)\)\s*;",text,re.S)
    if not match: raise RuntimeError(f"cannot parse field {path}")
    if vector: return [tuple(map(float,v)) for v in re.findall(r"\(([-+0-9.eE]+) ([-+0-9.eE]+) ([-+0-9.eE]+)\)",match.group(1))]
    return [float(v) for v in match.group(1).split()]

def labels(path):
    text=path.read_text(); match=re.search(r"\n\s*\d+\s*\((.*?)\)",text,re.S)
    return list(map(int,match.group(1).split())) if match else []

def percentile(data,p):
    ordered=sorted(data); x=(len(ordered)-1)*p; lo=int(x); hi=min(lo+1,len(ordered)-1)
    return ordered[lo]+(ordered[hi]-ordered[lo])*(x-lo)

def nearest_master_segment(y,z,contract):
    w=contract["dimensions"]["width"]; h=contract["dimensions"]["height"]; r=contract["upperFillets"]["radius"]
    # Exact Q1 chord junctions used by the generator: right wall, two right
    # chords, roof, two left chords, left wall.
    s=math.sqrt(0.5); chain=[(w,0),(w,h-r),(w-r+r*s,h-r+r*s),(w-r,h),(r,h),(r-r*s,h-r+r*s),(0,h-r),(0,0)]
    best=None
    for index,(a,b) in enumerate(zip(chain,chain[1:])):
        ay,az=a; by,bz=b; dy,dz=by-ay,bz-az; den=dy*dy+dz*dz
        t=max(0.0,min(1.0,((y-ay)*dy+(z-az)*dz)/den)); q=(ay+t*dy,az+t*dz)
        d=math.hypot(y-q[0],z-q[1])
        candidate=(d,index,q,t)
        if best is None or candidate<best: best=candidate
    return best

def main():
    parser=argparse.ArgumentParser(); parser.add_argument("candidate",type=Path); parser.add_argument("--contract",type=Path,default=Path("geometry/interface-geometry-v1.json")); parser.add_argument("--csv",type=Path,required=True); parser.add_argument("--json",type=Path,required=True)
    args=parser.parse_args(); mesh=args.candidate/"convert/constant/polyMesh"; fields=args.candidate/"convert/constant"
    det=values(fields/"cellDeterminant"); vol=values(fields/"cellVolume"); centres=values(fields/"C",True)
    owner=labels(mesh/"owner"); neighbour=labels(mesh/"neighbour"); adjacent=[set() for _ in det]
    for face,b in enumerate(neighbour):
        a=owner[face]; adjacent[a].add(b); adjacent[b].add(a)
    contract=json.loads(args.contract.read_text()); axial=contract["discretization"]["fluidAxialCells"]
    manifest=json.loads((args.candidate/"generation-manifest.json").read_text()); axial=manifest["effectiveDiscretization"]["fluidAxialCells"]
    cross=len(det)//axial; length=contract["dimensions"]["length"]
    rows=[]
    for cell,d in enumerate(det):
        if d>=1e-3: continue
        x,y,z=centres[cell]; distance,segment,q,t=nearest_master_segment(y,z,contract)
        axial_index=min(axial-1,max(0,int(x/length*axial)))
        rows.append({"cellId":cell,"centerX":x,"centerY":y,"centerZ":z,"volume":vol[cell],"determinant":d,
          "generatorBlock":"single-transfinite-extrusion","axialIndex":axial_index,"crossSectionIndex":cell//axial,
          "distanceToFilletedInterface":distance,"distanceToInterface":distance,"masterFaceSegment":segment,
          "axialCoordinate":x,"transverseLayer":"interface-adjacent" if distance<2e-6 else "interior",
          "neighbors":" ".join(map(str,sorted(adjacent[cell]))),"closestY":q[0],"closestZ":q[1],"segmentCoordinate":t})
    args.csv.parent.mkdir(parents=True,exist_ok=True)
    with args.csv.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]) if rows else ["cellId"],lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
    summary={"cellCount":len(det),"badCellCount":len(rows),"minimumDeterminant":min(det),"maximumDeterminant":max(det),
      "determinantPercentiles":{"p01":percentile(det,.01),"p05":percentile(det,.05),"p50":percentile(det,.5),"p95":percentile(det,.95),"p99":percentile(det,.99)},
      "cumulativeHistogram":{format(t,".0e"):sum(v<t for v in det) for t in THRESHOLDS},
      "axialIndicesWithBadCells":sorted({r["axialIndex"] for r in rows}),"masterSegmentsWithBadCells":{str(s):sum(r["masterFaceSegment"]==s for r in rows) for s in range(7)},
      "interfaceAdjacentBadCells":sum(r["transverseLayer"]=="interface-adjacent" for r in rows),
      "neighborVolumeRatio":{"maximum":max(max(vol[a],vol[b])/min(vol[a],vol[b]) for a,neighbors in enumerate(adjacent) for b in neighbors if a<b),
        "p99":percentile([max(vol[a],vol[b])/min(vol[a],vol[b]) for a,neighbors in enumerate(adjacent) for b in neighbors if a<b],.99)}}
    args.json.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n"); print(json.dumps(summary,sort_keys=True))

if __name__=="__main__": main()
