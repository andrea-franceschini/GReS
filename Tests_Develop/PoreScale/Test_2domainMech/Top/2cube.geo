Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

lc = 0.25;
// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:

Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};

Point(5) = {0, 0, 2, lc};
Point(6) = {1, 0, 2, lc};
Point(7) = {1, 1, 2, lc};
Point(8) = {0, 1, 2, lc};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Curve Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1};
Physical Surface("bot1",1)={1};

Plane Surface(2) = {2};
Physical Surface("bot2",2)={2};


Transfinite Surface {1};
Recombine Surface {1};
Transfinite Surface {2};
Recombine Surface {2};

Extrude {0, 0,0.5} { Surface{1}; Layers{4}; Recombine;}
Extrude {0, 0,0.5} { Surface{2}; Layers{4}; Recombine;}

Physical Volume("v1", 1) = {1};
Physical Volume("v2", 2) = {2};


Mesh 3;
Save "Test2blocks.msh";

