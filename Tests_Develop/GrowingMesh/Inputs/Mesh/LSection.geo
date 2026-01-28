Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
// Mesh.Format = 16; // vtk output format

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:
lc = 1e-2;
fileName = "domainL.msh";
// fileName = "domainL.vtk";

// Geometry
//  XZ view
//               7------8
//               |      |
//  Z(2)         |      |
//               |      |
//       4-------5------6
//       |       |      |
//  Z(1) |       |      |
//       1-------2------3
//          X(1)    X(2)

// Domain dimension
CompX(1) = 3;
CompX(2) = 10;
CompY    = 4;
CompZ(1) = 1;
CompZ(2) = 5;

ElmX[1] = 6;
ElmX[2] = 20;
ElmY    = 8;
ElmZ[1] = 2;
ElmZ[2] = 10;

// Domain creation
// Points
px=0;
py=0;
pz=0;
Point(1) = {px,py,pz,lc};
px=CompX(1);
Point(2) = {px,py,pz,lc};
px=CompX(1)+CompX(2);
Point(3) = {px,py,pz,lc};


px=0;
py=CompY;
Point(4) = {px,py,pz,lc};
px=CompX(1);
Point(5) = {px,py,pz,lc};
px=CompX(1)+CompX(2);
Point(6) = {px,py,pz,lc};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};

Line(3) = {4, 5};
Line(4) = {5, 6};

Line(5) = {1, 4};
Line(6) = {2, 5};
Line(7) = {3, 6};

Transfinite Line{1,3} = ElmX(1)+1;
Transfinite Line{2,4} = ElmX(2)+1;

Transfinite Line{5,6,7} = ElmY+1;

// Surfaces
Curve Loop(1) = {1, 6, -3, -5};
Curve Loop(2) = {2, 7, -4, -6};
For s In {1:2}
    Plane Surface(s) = {s};
    Transfinite Surface {s}; // structured grid
    Recombine Surface {s}; // using hexahedra
EndFor

Extrude {0, 0, CompZ(1)} { Surface{1:2}; Layers{ElmZ(1)}; Recombine;}
Extrude {0, 0, CompZ(2)} { Surface{29}; Layers{ElmZ(2)}; Recombine;}

Physical Volume("domain", 1) = {1:3};
Physical Surface("latX0",1) = {28,72};
Physical Surface("latXM",2) = {42};
Physical Surface("latY0",3) = {16,38,60};
Physical Surface("latYM",4) = {24,46,68};
Physical Surface("latZ0",5) = {1:2};
Physical Surface("latZM",6) = {51,64,73};

Physical Surface("FaceA",7) = {51};
Physical Surface("FaceB",8) = {64};

Mesh 3;
Save Sprintf(fileName);
