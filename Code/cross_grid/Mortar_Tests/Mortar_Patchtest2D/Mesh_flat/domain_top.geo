// Gmsh project created on Fri Mar 15 13:58:35 2024

Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

lc = 0.1;
// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:

Point(1) = {0, 1, 0, lc};
Point(2) = {1, 1, 0, lc};
Point(3) = {1, 2, 0, lc};
Point(4) = {0, 2, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};


Plane Surface(1) = {1};
Transfinite Surface {1};

Physical Curve("Interface_top2bot",1) = {1};
Physical Curve("Load_edge",2) = {3};
Physical Curve("Lateral_fixed",3) = {2,4};
Physical Surface("Domain_1",1)={1};

Mesh 2;
Save "TopBlock_tetra.msh";

