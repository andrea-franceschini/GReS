Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:
lc = 1e-2;

// Domain dimension
Dim[1] = 5;
Dim[2] = 10;
Dim[3] = 10;

Elms[1] = 10; // 2  16
Elms[2] = 40; // 2  16
Elms[3] = 40; // 4  32

// Fault partition
PlaEsp = 1;
PlaAlp = (90./180.)*Pi; // (30./90.)*Pi Pi/2
PlaGam = Pi/2.-PlaAlp;

compY = Cos(PlaAlp)*Dim[2];
// Domain creation
// Points
posx[1] = Dim[1];
posy[1] = Dim[2];
posz[1] = Dim[3];
posy[2] = posy[1]+((1-PlaPer)/2)*compY;

Point(1) = {0,      0,       0, lc};
Point(2) = {0,posy[1],       0, lc};
Point(3) = {0,posy[2], posz[1], lc};
Point(4) = {0,      0, posz[1], lc};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line{1,3} = Elms[2]+1;
Transfinite Line{2,4} = Elms[3]+1;

Curve Loop(1) = { 1, 2, 3, 4};

For s In {1:1}
    Plane Surface(s) = {s};
    Transfinite Surface {s}; // structured grid
    Recombine Surface {s}; // using hexahedra
EndFor

Extrude {Dim[1], 0, 0} { Surface{1}; Layers{Elms[1]}; Recombine;}

Physical Volume("domain", 1) = {1};
Physical Surface("latX0",1) = {1};
Physical Surface("latXM",2) = {26};
Physical Surface("latY0",3) = {25};
Physical Surface("latYM",4) = {17};
Physical Surface("latZ0",5) = {13};
Physical Surface("latZM",6) = {21};

Mesh 3;
Save "Fault.msh";