Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:
lc = 1e-2;

// Domain dimension
H = 10.;
h = 5.;
L = 10.;
l = 5.;

DimX[1] = 0.;
DimX[2] = (H-h)/2.;
DimX[3] = DimX[2]+h;
DimX[4] = H;

DimY[1] = 0.;
DimY[2] = (L-l)/2.;
DimY[3] = DimY[2]+l;
DimY[4] = L;

DimZ[1] = 2.;
DimZ[2] = 4+DimY[2];
DimZ[3] = 2+DimZ[3];

ElmsX[1] = 8; // 2  16
ElmsX[2] = 16; // 2  16
ElmsX[3] = 8; // 4  32

ElmsY[1] = 8; // 2  16
ElmsY[2] = 16; // 2  16
ElmsY[3] = 8; // 4  32

ElmsZ[1] = 4; // 2  16
ElmsZ[2] = 8; // 2  16
ElmsZ[3] = 4; // 4  32

// Domain creation
// Points
For s In {1:4}
    Point(4*(s-1)+1) = {DimX[1],DimY[s],DimZ[1], lc};
    Point(4*(s-1)+2) = {DimX[2],DimY[s],DimZ[1], lc};
    Point(4*(s-1)+3) = {DimX[3],DimY[s],DimZ[1], lc};
    Point(4*(s-1)+4) = {DimX[4],DimY[s],DimZ[1], lc};
EndFor

// Lines
For s In {1:4}
    Line(3*(s-1)+1) = {4*(s-1)+1, 4*(s-1)+2};
    Line(3*(s-1)+2) = {4*(s-1)+2, 4*(s-1)+3};
    Line(3*(s-1)+3) = {4*(s-1)+3, 4*(s-1)+4};
EndFor

For s In {1:4}
    Line(3*(s-1)+13) = {  s, s+4};
    Line(3*(s-1)+14) = {s+4, s+8};
    Line(3*(s-1)+15) = {s+8,s+12};
EndFor

Curve Loop(1) = { 1, 16,-4,-13};
Curve Loop(2) = { 2, 19,-5,-16};
Curve Loop(3) = { 3, 22,-6,-19};

Curve Loop(4) = { 4, 17,-7,-14};
Curve Loop(5) = { 5, 20,-8,-17};
Curve Loop(6) = { 6, 23,-9,-20};

Curve Loop(7) = { 7, 18,-10,-15};
Curve Loop(8) = { 8, 21,-11,-18};
Curve Loop(9) = { 9, 24,-12,-21};

// Mesh subdivision.
Transfinite Line{1,4,7,10} = ElmsX[1]+1;
Transfinite Line{2,5,8,11} = ElmsX[2]+1;
Transfinite Line{3,6,9,12} = ElmsX[3]+1;

Transfinite Line{13,16,19,22} = ElmsY[1]+1;
Transfinite Line{14,17,20,23} = ElmsY[2]+1;
Transfinite Line{15,18,21,24} = ElmsY[3]+1;

// Creating the planes.
For s In {1:9}
    Plane Surface(s) = {s};
    Transfinite Surface {s}; // structured grid
    Recombine Surface {s}; // using hexahedra
EndFor

// Extruding the surfaces.
Extrude {0, 0, DimZ(1)} { Surface{1:9}; Layers{ElmsZ[1]}; Recombine;}
Extrude {0, 0, DimZ(2)} { Surface{46,68,90,112,134,156,178,200,222}; Layers{ElmsZ[2]}; Recombine;}
Extrude {0, 0, DimZ(3)} { Surface{332}; Layers{ElmsZ[3]}; Recombine;}

// Select reigions.

Physical Volume("Soil", 1) = {1:13,15,16,17,18};
Physical Volume("Build", 2) = {14,19};
Physical Surface("latZ0",5) = {1:9};
/*
Physical Surface("latX0",1) = { 1:9 };
Physical Surface("latXM",2) = {47,69,91,113,135,157,179,201,223};
Physical Surface("latY0",3) = {46,68,90};
Physical Surface("latYM",4) = {170,192,214};
Physical Surface("latZM",6) = {86,152,218};

Mesh 3;
Save "Fault.msh";