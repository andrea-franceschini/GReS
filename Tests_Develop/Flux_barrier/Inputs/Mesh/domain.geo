Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
// Mesh.Format = 16; // vtk output format

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:
lc = 1e-2;
fileName = "Fault.msh";
// fileName = "Fault.vtk";

// Domain dimension
Dim[1] = 1;
Dim[2] = 1;
Dim[3] = 1;

Elms[1] = 8;
Elms[2] = 8;
Elms[3] = 8;

// Fault partition
PlaDim = 0.5;
PlaEsp = 0.1;
PlaAlp = (30./90.)*Pi;
PlaGam = Pi/2.-PlaAlp;

ElmSF[1] = 4;
ElmSF[2] = 4;
ElmSF[3] = Elms[3];

compY = Cos(PlaAlp)*Dim[2];
// Domain creation
// Points
posx[1] = Dim[1];
posy[1] = Dim[2];
posz[1] = Dim[3];
posy[2] = posy[1]+compY;

Point(1) = {0,      0,       0, lc};
Point(2) = {0,posy[1],       0, lc};
Point(3) = {0,posy[2], posz[1], lc};
Point(4) = {0,      0, posz[1], lc};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line{1,3} = Elms[2]+1;
Transfinite Line{2,4} = Elms[3]+1;

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Surface {1}; // structured grid
Recombine Surface {1}; // using hexahedra

// Fault creation
posy[3] = posy[1]+PlaEsp;
posy[4] = posy[2]+PlaEsp;
pt = newp;
Point(pt+1) = {0,posy[3],       0, lc};
Point(pt+2) = {0,posy[4], posz[1], lc};

// Lines
nl = newl;
Line(nl+1) = {2,pt+1};
Line(nl+2) = {pt+1,pt+2};
Line(nl+3) = {pt+2,3};

Transfinite Line{nl+1,nl+3} = ElmSF[2]+1;
Transfinite Line{nl+2} = ElmSF[3]+1;

Curve Loop(2) = {nl+1,nl+2,nl+3,-2};
Plane Surface(2) = {2};

Transfinite Surface {2}; // structured grid
Recombine Surface {2}; // using hexahedra

// Domain creation
posy[5] = posy[3]+Dim[2]+compY;
posy[6] = posy[4]+Dim[2];
pt = newp;
Point(pt+1) = {0,posy[5],       0, lc};
Point(pt+2) = {0,posy[6], posz[1], lc};

// Lines
nl = newl;
Line(nl+1) = {6,pt+1};
Line(nl+2) = {pt+1,pt+2};
Line(nl+3) = {pt+2,7};

Transfinite Line{nl+1,nl+3} = Elms[2]+1;
Transfinite Line{nl+2} = Elms[3]+1;

Curve Loop(3) = {nl+1,nl+2,nl+3,-7};
Plane Surface(3) = {3};

Transfinite Surface {3}; // structured grid
Recombine Surface {3}; // using hexahedra

Extrude {Dim[1], 0, 0} { Surface{1:3}; Layers{Elms[1]}; Recombine;}

Physical Volume("domain", 1) = {1,3};
Physical Volume("core", 2) = {2};
Physical Surface("latX0",1) = { 1, 2, 3};
Physical Surface("latXM",2) = {34,56,78};
Physical Surface("latY0",3) = {33};
Physical Surface("latYM",4) = {69};
Physical Surface("latZ0",5) = {21,43,65};
Physical Surface("latZM",6) = {29,51,73};

Mesh 3;
Save Sprintf(fileName);
