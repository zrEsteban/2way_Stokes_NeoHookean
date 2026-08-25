SetFactory("OpenCASCADE");
Mesh.RecombineAll = 1;
Rectangle(1) = {0, 0, 0, 1e-3, 2e-4};
out[] = Extrude {0, 0, 1e-4} { Surface{1}; Layers{2}; Recombine; };
// Extrusion returns top, volume, then four lateral surfaces.
Physical Surface(1) = {out[2], out[3], out[4], out[5]};
Physical Surface(2) = {1};
Physical Surface(3) = {out[0]};
Physical Volume(10) = {out[1]};
