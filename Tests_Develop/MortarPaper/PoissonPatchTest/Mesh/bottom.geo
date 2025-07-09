Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

nX = 4;
nY = 4;
nZ = 4;
// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:

X = 1;
Y = 1;
Z = 1;
xmin = 0;
ymin = 0;
zmin = 0;

Point(1) = {xmin, ymin, zmin};
Point(2) = {xmin+X, ymin, zmin};
Point(3) = {xmin+X, ymin+Y, zmin};
Point(4) = {xmin, ymin+Y, zmin};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line{1} = nX+1;
Transfinite Line{2} = nY+1;
Transfinite Line{3} = nX+1;
Transfinite Line{4} = nY+1; 

Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Transfinite Surface {1}; // structured grid
Recombine Surface {1}; // using hexahedra

// extruding mesh along the z direction
Extrude {0, 0, Z} { Surface{1}; Layers{nZ}; Recombine;}
Physical Volume("cubeBot",1) = {1};
Physical Surface("lat",1) = {1,13,17,21,25};
Physical Surface("interf",2) = {26};

// Mesh second-order
Mesh.ElementOrder = 1;
//Mesh.SecondOrderIncomplete = 0;  // <== ENSURES FULL HEXA-27

Mesh 3;
Save "bottom.msh";