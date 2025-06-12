Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:
lc = 1e-2;

// Models dimension
x_len = 5; // 0.1 / 1
y_len = 5; // 0.1 / 1
z_len = 5; // 0.1 / 10

// Mesh partition
x_elm = 20; //  1 /  4
y_elm = 20; //  1 /  4
z_elm = 20; // 30 / 40

// Model creation
// Points
Point(4) = {-x_len/2.,-y_len/2., 0, lc};
Point(1) = { x_len/2.,-y_len/2., 0, lc};
Point(2) = { x_len/2., y_len/2., 0, lc};
Point(3) = {-x_len/2., y_len/2., 0, lc};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line{1} = y_elm+1;
Transfinite Line{2} = x_elm+1;
Transfinite Line{3} = y_elm+1;
Transfinite Line{4} = x_elm+1;

Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Transfinite Surface {1}; // structured grid
Recombine Surface {1}; // using hexahedra

// extruding mesh along the z direction
Extrude {0, 0, z_len} { Surface{1}; Layers{z_elm}; Recombine;}

Physical Volume("column", 1) = {1};
Physical Surface("top",2) = {26};
Physical Surface("bot",1) = {1};
Physical Surface("latX0",5) = {21};
Physical Surface("latY0",6) = {25};
Physical Surface("latXM",3) = {13};
Physical Surface("latYM",4) = {17};

Mesh 3;
Save "Column.msh";
