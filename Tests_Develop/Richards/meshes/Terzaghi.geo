Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

lc = 0.125;
// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:

Point(1) = {-0.5, -0.5, 0, lc};
Point(2) = {0.5,-0.5, 0, lc};
Point(3) = {0.5, 0.5, 0, lc};
Point(4) = {-0.5, 0.5, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};
Physical Surface("dir_bottom",1)={1};




Transfinite Surface {1};
Recombine Surface {1};

Extrude {0, 0, 10} { Surface{1}; Layers{80}; Recombine;}
Physical Volume("Sub1", 1) = {1};
Mesh 3;
Save "Richards_refined.msh";

