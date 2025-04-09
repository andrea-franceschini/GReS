// Gmsh project created on Fri Mar 15 13:58:35 2024

Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

//lc = 0.04;
// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:
nX = 16;
nY = 16;

Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 0.5, 0};
Point(4) = {0, 0.5, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line{1} = nX;
Transfinite Line{2} = nY;
Transfinite Line{3} = nX;
Transfinite Line{4} = nY;
Transfinite Surface {1};
Recombine Surface {1};

Physical Curve("Interface_bot2top",1) = {3};
Physical Curve("Fixed_bottom",2) = {1};
Physical Curve("Lateral_fixed",3) = {2};
Physical Surface("Domain_2",1)={1};

//Mesh 2;
//Save "bottom.msh";

