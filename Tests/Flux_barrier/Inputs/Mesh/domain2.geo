Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:
lc = 1e-2;

// Domain dimension
Dim[1] = 10;
Dim[2] = 10;
Dim[3] = 10;

Elms[1] = 16; // 2  16
Elms[2] = 16; // 2  16
Elms[3] = 32; // 4  32

// Fault partition
PlaPer = 0.5;  // percent
PlaEsp = 1;
PlaAlp = Pi/2; // (30./90.)*Pi Pi/2
PlaGam = Pi/2.-PlaAlp;

ElmSF[1] = 16; // 2   16
ElmSF[2] = 8; // 1    8
ElmSF[3] = Elms[3];

compY = Cos(PlaAlp)*Dim[2];
// Domain creation
// Points
posx[1] = Dim[1];
posy[1] = Dim[2];
posz[1] = ((1-PlaPer)/2)*Dim[3];
posy[2] = posy[1]+((1-PlaPer)/2)*compY;

Point(1) = {0,      0,       0, lc};
Point(2) = {0,posy[1],       0, lc};
Point(3) = {0,posy[2], posz[1], lc};
Point(4) = {0,      0, posz[1], lc};

posy[3] = posy[2]+PlaPer*compY;
posz[2] = posz[1]+PlaPer*Dim[3];
Point(5) = {0,posy[3], posz[2], lc};
Point(6) = {0,      0, posz[2], lc};

posy[4] = posy[1]+compY;
posz[3] = posz[2]+((1-PlaPer)/2)*Dim[3];
Point(7) = {0,posy[4], posz[3], lc};
Point(8) = {0,      0, posz[3], lc};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(6) = {3, 5};
Line(7) = {5, 6};
Line(8) = {6, 4};

Line( 9) = {5, 7};
Line(10) = {7, 8};
Line(11) = {8, 6};

Transfinite Line{1,3,7,10} = Elms[2]+1;
Transfinite Line{6,8} = Round(PlaPer*Elms[3])+1;
Transfinite Line{2,4,9,11} = Round(((1-PlaPer)/2)*Elms[3])+1;

Curve Loop(1) = { 1, 2, 3, 4};
Curve Loop(2) = {-3, 6, 7, 8};
Curve Loop(3) = {-7, 9,10,11};

For s In {1:3}
    Plane Surface(s) = {s};
    Transfinite Surface {s}; // structured grid
    Recombine Surface {s}; // using hexahedra
EndFor

// Fault creation
posy[5] = posy[1]+PlaEsp;
posy[6] = posy[2]+PlaEsp;
Point( 9) = {0,posy[5],       0, lc};
Point(10) = {0,posy[6], posz[1], lc};

posy[7] = posy[3]+PlaEsp;
Point(11) = {0,posy[7], posz[2], lc};

posy[8] = posy[4]+PlaEsp;
Point(12) = {0,posy[8], posz[3], lc};

// Lines
Line(12) = { 2, 9};
Line(13) = { 9,10};
Line(14) = {10,3};

Line(15) = {10,11};
Line(16) = {11, 5};

Line(17) = {11,12};
Line(18) = {12,7};

Transfinite Line{12,14,16,18} = ElmSF[2]+1;
Transfinite Line{15} = Round(PlaPer*ElmSF[3])+1;
Transfinite Line{13,17} = Round(((1-PlaPer)/2)*ElmSF[3])+1;

Curve Loop(4) = { 12,13,14,-2};
Curve Loop(5) = {-14,15,16,-6};
Curve Loop(6) = {-16,17,18,-9};

For s In {4:6}
    Plane Surface(s) = {s};
    Transfinite Surface {s}; // structured grid
    Recombine Surface {s}; // using hexahedra
EndFor

// Domain creation
posy[ 9] = posy[5]+Dim[2]+compY;
posy[10] = posy[6]+Dim[2]+((1-PlaPer)/2+PlaPer)*compY;
Point(13) = {0,posy[ 9],      0,lc};
Point(14) = {0,posy[10],posz[1],lc};

posy[11] = posy[7]+Dim[2]+((1-PlaPer)/2)*compY;
Point(15) = {0,posy[11],posz[2],lc};

posy[12] = posy[8]+Dim[2];
Point(16) = {0,posy[12],posz[3],lc};

// Lines
Line(19) = { 9,13};
Line(20) = {13,14};
Line(21) = {14,10};

Line(22) = {14,15};
Line(23) = {15,11};

Line(24) = {15,16};
Line(25) = {16,12};

Transfinite Line{19,21,23,25} = Elms[2]+1;
Transfinite Line{22} = Round(PlaPer*Elms[3])+1;
Transfinite Line{20,24} = Round(((1-PlaPer)/2)*Elms[3])+1;

Curve Loop(7) = { 19,20,21,-13};
Curve Loop(8) = {-21,22,23,-15};
Curve Loop(9) = {-23,24,25,-17};

For s In {7:9}
    Plane Surface(s) = {s};
    Transfinite Surface {s}; // structured grid
    Recombine Surface {s}; // using hexahedra
EndFor

Extrude {Dim[1], 0, 0} { Surface{1:9}; Layers{Elms[1]}; Recombine;}

Physical Volume("domain", 1) = {1:4,6:9};
Physical Volume("core", 2) = {5};
Physical Surface("latX0",1) = { 1:9 };
Physical Surface("latXM",2) = {47,69,91,113,135,157,179,201,223};
Physical Surface("latY0",3) = {46,68,90};
Physical Surface("latYM",4) = {170,192,214};
Physical Surface("latZ0",5) = {34,100,166};
Physical Surface("latZM",6) = {86,152,218};

Mesh 3;
Save "Fault.msh";