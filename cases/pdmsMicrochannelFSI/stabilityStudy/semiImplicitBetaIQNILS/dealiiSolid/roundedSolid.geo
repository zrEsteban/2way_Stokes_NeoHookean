// Matching PDMS solid around a channel with rounded upper corners (SI units).
SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 2.2;
Mesh.RecombineAll = 1;
Mesh.Algorithm = 6;
Mesh.RecombinationAlgorithm = 1;
If (!Exists(allHex)) allHex = 0; EndIf
If (allHex)
  // deal.II 9.3 GridIn cannot ingest extruded prisms. Subdivide any residual
  // triangles into quadrilaterals before extrusion so the volume is all-hex.
  Mesh.SubdivisionAlgorithm = 1;
  Mesh.RecombinationAlgorithm = 2;
EndIf

If (!Exists(L)) L = 0.00492; EndIf
If (!Exists(w)) w = 0.000492; EndIf
If (!Exists(h)) h = 0.000033; EndIf
If (!Exists(r)) r = 0.00000075; EndIf
If (!Exists(t)) t = 0.0005; EndIf
If (!Exists(e)) e = 0.0005; EndIf
If (!Exists(nx)) nx = 75; EndIf
If (!Exists(lcNear)) lcNear = 0.00000075; EndIf
If (!Exists(lcFar)) lcFar = 0.00006; EndIf

Point(1) = {0, -t, 0, lcFar};
Point(2) = {0, -t, h+e, lcFar};
Point(3) = {0, w+t, h+e, lcFar};
Point(4) = {0, w+t, 0, lcFar};
Point(5) = {0, w, 0, lcFar};
Point(6) = {0, w, h-r, lcNear};
Point(7) = {0, w-r, h, lcNear};
Point(8) = {0, r, h, lcNear};
Point(9) = {0, 0, h-r, lcNear};
Point(10) = {0, 0, 0, lcFar};
Point(11) = {0, w-r, h-r, lcNear};
Point(12) = {0, r, h-r, lcNear};

Line(1) = {1, 2};             // outer left
Line(2) = {2, 3};             // outer roof
Line(3) = {3, 4};             // outer right
Line(4) = {4, 5};             // right base
Line(5) = {5, 6};             // right interface
Circle(6) = {6, 11, 7};       // right upper fillet
Line(7) = {7, 8};             // roof interface
Circle(8) = {8, 12, 9};       // left upper fillet
Line(9) = {9, 10};            // left interface
Line(10) = {10, 1};           // left base
Curve Loop(1) = {1,2,3,4,5,6,7,8,9,10};
Plane Surface(1) = {1};

Field[1] = Distance;
Field[1].CurvesList = {6, 8};
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lcNear;
Field[2].SizeMax = lcFar;
Field[2].DistMin = 2*r;
Field[2].DistMax = 4*r;
Background Field = 2;

Recombine Surface {1};
out[] = Extrude {L, 0, 0} { Surface{1}; Layers{nx}; Recombine; };

// Stable numeric tags are shared by OpenFOAM and the deal.II participant.
Physical Volume("solid", 10) = {out[1]};
Physical Surface("ends", 1) = {1, out[0]};
Physical Surface("outer", 2) = {out[2], out[3], out[4]};
Physical Surface("base", 3) = {out[5], out[11]};
Physical Surface("interface", 4) = {out[6], out[7], out[8], out[9], out[10]};
