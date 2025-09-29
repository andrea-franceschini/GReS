Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:
lc = 1e-2;

// Domain dimension
Dim[1] = 1;
Dim[2] = 1;
Dim[3] = 1;

Elms[1] = 2;
Elms[2] = 2;
Elms[3] = 2;

// Domain creation
// Points
posx[1] = Dim[1];
posy[1] = Dim[2];
posz[1] = Dim[3];

Point(1) = {     0,      0, 0, lc};
Point(2) = {Dim[1],      0, 0, lc};
Point(3) = {Dim[1], Dim[2], 0, lc};
Point(4) = {     0, Dim[2], 0, lc};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line{1,3} = Elms[1]+1;
Transfinite Line{2,4} = Elms[2]+1;

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Surface {1}; // structured grid
Recombine Surface {1}; // using hexahedra

Extrude {0, 0, Dim[3]} { Surface{1}; Layers{Elms[3]}; Recombine;}

Physical Volume("domain", 1) = {1};
Physical Surface("latX0",1) = {25};
Physical Surface("latXM",2) = {17};
Physical Surface("latY0",3) = {13};
Physical Surface("latYM",4) = {21};
Physical Surface("latZ0",5) = {1};
Physical Surface("latZM",6) = {26};

Mesh 3;
Save "Fault.msh";