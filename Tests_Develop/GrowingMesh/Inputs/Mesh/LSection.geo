Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
// Mesh.Format = 16; // vtk output format

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:
lc = 1e-2;
// fileName = "domainL.msh";
fileName = "domainL.vtk";

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
CompX(1) = 10;
CompX(2) = 1;
CompY    = 1;
CompZ(1) = 1;
CompZ(2) = 5;

ElmX[1] = 2;
ElmX[2] = 1;
ElmY    = 1;
ElmZ[1] = 1;
ElmZ[2] = 1;

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

pz=CompZ(1);
px=0;
Point(4) = {px,py,pz,lc};
px=CompX(1);
Point(5) = {px,py,pz,lc};
px=CompX(1)+CompX(2);
Point(6) = {px,py,pz,lc};

pz=CompZ(1)+CompZ(2);
px=CompX(1);
Point(7) = {px,py,pz,lc};
px=CompX(1)+CompX(2);
Point(8) = {px,py,pz,lc};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};

Line(3) = {4, 5};
Line(4) = {5, 6};

Line(5) = {7, 8};

Line(6) = {1, 4};
Line(7) = {2, 5};
Line(8) = {3, 6};

Line( 9) = {5, 7};
Line(10) = {6, 8};

Transfinite Line{1,3} = ElmX(1)+1;
Transfinite Line{2,4,5} = ElmX(2)+1;

Transfinite Line{6,7,8} = ElmZ(1)+1;
Transfinite Line{9,10} = ElmZ(2)+1;

// Surfaces
Curve Loop(1) = {1, 7, -3, -6};
Curve Loop(2) = {2, 8, -4, -7};
Curve Loop(3) = {4,10, -5, -9};
For s In {1:3}
    Plane Surface(s) = {s};
    Transfinite Surface {s}; // structured grid
    Recombine Surface {s}; // using hexahedra
EndFor

Extrude {0, CompY, 0} { Surface{1:3}; Layers{ElmY}; Recombine;}

Physical Volume("domain", 1) = {1:3};
Physical Surface("latX0",1) = {31};
Physical Surface("latXM",2) = {45,67};
Physical Surface("latY0",3) = {1:3};
Physical Surface("latYM",4) = {32,54,76};
Physical Surface("latZ0",5) = {19,41};
Physical Surface("latZM",6) = {27,71,75};

Physical Surface("FaceA",7) = {27};
Physical Surface("FaceB",8) = {75};

Mesh 3;
Save Sprintf(fileName);
/*