Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:
lc = 1e-2;

// Models dimension
x_len = 200;
y_len = 10;
z_len = 200;

int_rg_prop = 0.2;
x_int = int_rg_prop*x_len;
z_int = int_rg_prop*z_len;

// Mesh partition
x_elm =  20;
y_elm =  10;
z_elm =  20;

// Model creation
// Points
Point(1) = {-x_int/2, 0,-z_int/2, lc};
Point(2) = { x_int/2, 0,-z_int/2, lc};
Point(3) = { x_int/2, 0, z_int/2, lc};
Point(4) = {-x_int/2, 0, z_int/2, lc};

Point(5) = {-x_len/2, 0,-z_len/2, lc};
Point(6) = { x_len/2, 0,-z_len/2, lc};
Point(7) = { x_len/2, 0, z_len/2, lc};
Point(8) = {-x_len/2, 0, z_len/2, lc};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Transfinite Line{1,3} = x_elm+1;
Transfinite Line{2,4} = z_elm+1;

Transfinite Line{5,7} = x_elm+1;
Transfinite Line{6,8} = z_elm+1;

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1};
Plane Surface(2) = {1,2};

//Transfinite Surface {1,2}; // structured grid
Transfinite Surface {2}; // structured grid
Recombine Surface {1,2}; // using hexahedra

/*
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
Save "Beam.msh";
