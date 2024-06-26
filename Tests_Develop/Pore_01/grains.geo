Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
SetFactory("OpenCASCADE");
lc = 0.02;
// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:

Point(1) = {0, 0, 0, lc};
Point(2) = {0.3, 0, 0, lc};
Point(3) = {0.3, 1, 0, lc};
Point(4) = {0, 1, 0, lc};
Point(5) = {-0.4,0.5,0,lc};
Line(1) = {1, 2};
Circle(2) = {2,5,3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};
//Physical Surface("bot",1)={1};
//Transfinite Surface {1};
Recombine Surface {1};
Extrude {0, 0,0.1} { Surface{1}; Layers{3}; Recombine;}
Physical Volume("g_left", 1) = {1};
Physical Surface("fixL",1)={25};
Physical Surface("mortarL",2)={17};

Point(17) = {0.7, 0, 0, lc};
Point(18) = {1, 0, 0, lc};
Point(19) = {1, 1, 0, lc};
Point(20) = {0.7, 1, 0, lc};
Point(21) = {1.4,0.5,0,lc};
Line(21) = {17, 18};
Line(22) = {18,19};
Line(23) = {19,20};
Circle(24) = {20,21,17};

Curve Loop(25) = {21, 22, 23, 24};

Plane Surface(27) = {25};
//Physical Surface("bot",1)={1};
//Transfinite Surface {1};
Recombine Surface {27};
Extrude {0,0,0.1} { Surface{27}; Layers{3}; Recombine;}
Physical Volume("g_right", 2) = {2};
Physical Surface("fixR",3)={40};
Physical Surface("mortarR",4)={48};
//
Box(3) = {0,0,0,1,1,0.1};
MeshSize{PointsOf{ Volume{3};}} = 0.005;
BooleanDifference{ Volume{3}; Delete;}{ Volume{1}; Delete;}
BooleanDifference{ Volume{3}; Delete;}{ Volume{2}; Delete;}
Physical Volume("fluid", 3) = {3};
Recombine Surface {33,34};
Mesh 3;
Save "FlowMatrix.msh";

