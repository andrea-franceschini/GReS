//Mesh.Format = 16;
Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

lc = 10;

// DISTANCES__________________________________________
// X - 590 m tot - building space 390m
H = 100. ;
dx = 30. ;
h = dx*13 + H*2 ;

// Y - 860 m tot - building space 660 m
V = 100. ;
dy = 30. ;
v = dy*22 + V*2 ;

// COORDINATES________________________________________

// X coordinates
X[1] = 0. ;
X[2] = h ;
X[3] = h ;
X[4] = 0.;

// Y coordinates
Y[1] = 0. ;
Y[2] = 0. ;
Y[3] = v ;
Y[4] = v ;

// Z coordinates
Z[0] = 0. ;


Z[1] = 15.48;
Z[2] = 3.33;
Z[3] = 2.22;
// Z[4] = 5.65;
// Z[5] = 2.14;
// Z[6] = 6.19;
// Z[7] = 3.57;
// Z[8] = 10.48;
// Z[9] = 1.90;

//MESH SUBDIVISION_________________________________

// like this each square is 29.5 x 30 m
//orthogonal to X
ElmsX = 20. ; //22
//orthogonal to Y
ElmsY = 30. ; //32

//orthogonal to Z
ElmsZ[1] = 2;  //31
ElmsZ[2] = 2;   //7
ElmsZ[3] = 2;   //4
// ElmsZ[4] = 11;
// ElmsZ[5] = 4;
// ElmsZ[6] = 12;
// ElmsZ[7] = 7;
// ElmsZ[8] = 21;
// ElmsZ[9] = 4;

// POINTS_________________________________________

Point(1) = {X[1], Y[1], Z[0], lc} ;
Point(2) = {X[2], Y[2], Z[0], lc} ;
Point(3) = {X[3], Y[3], Z[0], lc} ;
Point(4) = {X[4], Y[4], Z[0], lc} ;

// LINES__________________________________________

Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;

// TRANSFINITE LINES_____________________________

// Orthogonal to X
Transfinite Line{1,3} = ElmsX + 1 ; 

// Orthogonal to Y
Transfinite Line{2,4} = ElmsY + 1 ;

// CURVE LOOPS___________________________________

Curve Loop(1) = {1,2,3,4} ;
Plane Surface(1) = 1 ;

Transfinite Surface(1);
Recombine Surface(1);

// EXTRUSION_____________________________________

// Layer 1
Extrude{0,0,Z[1]} { Surface{1}; Layers{ElmsZ[1]}; Recombine;}

// Layer 2
Extrude{0,0,Z[2]} { Surface{26}; Layers{ElmsZ[2]}; Recombine;}

// Layer 3
Extrude{0,0,Z[3]} { Surface{48}; Layers{ElmsZ[3]}; Recombine;}

// PHYSICAL VOLUMES_____________________________

Physical Volume("Layer_1",1) = {1};
Physical Volume("Layer_2",2) = {2};
Physical Volume("Layer_3",3) = {3};

// PHYSICAL SURFACES____________________________

Physical Surface("X_0",1) = {25,47,69};
Physical Surface("X_M",2) = {17,39,61};
Physical Surface("Y_0",3) = {13,35,57};
Physical Surface("Y_M",4) = {21,43,65};
Physical Surface("Z_0",5) = {1}; 
Physical Surface("Z_M",6) = {70}; // TOP INTERFACE

Mesh 3;
Save "Bottom_part.msh" ;
//Save "Bottom_part.vtk" ;