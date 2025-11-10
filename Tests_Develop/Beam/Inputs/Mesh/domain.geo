Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
// Mesh.Format = 16; // vtk output format

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:
lc = 1e-2;
fileName = "Beam.msh";
// fileName = "Beam.vtk";

// Models dimension
x_len = 1.;
y_len = 0.1;
z_len = 0.1;

// Mesh partition
x_elm =  10;
y_elm =  1;
z_elm =  1;

// Model creation
// Points
Point(1) = {    0,     0, 0, lc};
Point(2) = {x_len,     0, 0, lc};
Point(3) = {x_len, y_len, 0, lc};
Point(4) = {    0, y_len, 0, lc};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line{1,3} = x_elm+1;
Transfinite Line{2,4} = y_elm+1;

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Surface {1}; // structured grid
Recombine Surface {1}; // using hexahedra

// extruding mesh along the z direction
Extrude {0, 0, z_len} { Surface{1}; Layers{z_elm}; Recombine;}

Physical Volume("column", 1) = {1};
Physical Surface("latX0",1) = {25};
Physical Surface("latXM",2) = {17};
Physical Surface("latY0",3) = {13};
Physical Surface("latYM",4) = {21};
Physical Surface("latZ0",5) = {1};
Physical Surface("latZM",6) = {26};

Mesh 3;
Save Sprintf(fileName);