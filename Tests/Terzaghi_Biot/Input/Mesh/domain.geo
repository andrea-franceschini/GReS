Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
// Mesh.Format = 16; // VTK output format

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:
lc = 1e-2;
meshFile = "Column_hexa.msh";
//meshFile = "Column.msh";

Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line{1} = 3;
Transfinite Line{2} = 3;
Transfinite Line{3} = 3;
Transfinite Line{4} = 3;

Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Transfinite Surface {1}; // structured grid
Recombine Surface {1}; // using hexahedra

// extruding mesh along the z direction
Extrude {0, 0, 10} { Surface{1}; Layers{20}; Recombine;}

Physical Volume("column", 1) = {1};
Physical Surface("top",2) = {26};
Physical Surface("bot",1) = {1};
Physical Surface("latX",4) = {17,25};
Physical Surface("latY",3) = {13,21};

Mesh 3;
Save Sprintf(meshFile);