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

// Models dimension (double the size)
x_len = 200;
y_len = 200;
z_len = 10;

radius = 15;
x_int = 10;
y_int = 10;

// Mesh partition
x_elm = 20;
y_elm = 20;
z_elm = 10;
r_elm = 10;

// Model creation
// Points
Point( 1) = {-x_int,-y_int, 0, lc};
Point( 2) = { x_int,-y_int, 0, lc};
Point( 3) = { x_int, y_int, 0, lc};
Point( 4) = {-x_int, y_int, 0, lc};

Point( 5) = {-x_int-radius,-y_int,        0, lc};
Point( 6) = {       -x_int,-y_int-radius, 0, lc};
Point( 7) = {        x_int,-y_int-radius, 0, lc};
Point( 8) = { x_int+radius,-y_int,        0, lc};
Point( 9) = { x_int+radius, y_int,        0, lc};
Point(10) = {        x_int, y_int+radius, 0, lc};
Point(11) = {       -x_int, y_int+radius, 0, lc};
Point(12) = {-x_int-radius, y_int,        0, lc};

Point(13) = {-x_len,-y_len, 0, lc};
Point(14) = { x_len,-y_len, 0, lc};
Point(15) = { x_len, y_len, 0, lc};
Point(16) = {-x_len, y_len, 0, lc};

// Lines
Circle(1) = {5,1,6};
Line(2) = {6, 7};
Circle(3) = {7,2,8};
Line(4) = {8, 9};
Circle(5) = {9,3,10};
Line(6) = {10,11};
Circle(7) = {11,4,12};
Line(8) = {12, 5};

Line( 9) = {13,14};
Line(10) = {14,15};
Line(11) = {15,16};
Line(12) = {16,13};

Transfinite Line{4,8} = Round(y_elm/2)+1;
Transfinite Line{2,6} = Round(x_elm/2)+1;
Transfinite Line{1,3,5,7} = r_elm+1;

Transfinite Line{10,12} = y_elm+1;
Transfinite Line{ 9,11} = x_elm+1;

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Curve Loop(2) = {9,10,11,12};
Plane Surface(1) = {1};
Plane Surface(2) = {1,2};

//Transfinite Surface {1,2}; // structured grid
//Transfinite Surface {1}; // structured grid
Recombine Surface {1,2}; // using hexahedra

// extruding mesh along the z direction
Extrude {0,0,z_len} { Surface{1,2}; Layers{z_elm}; Recombine;}

Physical Volume("volume", 1) = {1,2};

Physical Surface("latX1",1) = {103};
Physical Surface("latX2",2) = {111};
Physical Surface("latY1",3) = {115};
Physical Surface("latY2",4) = {107};
Physical Surface("latZ1",5) = {1,2};
Physical Surface("latZOUT",6) = {54};
Physical Surface("latZIN",7) = {116};

/*
Mesh 3;
Save Sprintf(fileName);
*/
