Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

lc = 0.1;

// COORDINATES________________________________________

// X coordinates
X[1] = 0. ;
X[2] = 100. ;
X[3] = 100. ;
X[4] = 0.;

// Y coordinates
Y[1] = 0. ;
Y[2] = 0. ;
Y[3] = 150. ;
Y[4] = 150. ;

// Z coordinates
Z[0] = 30. ;

Z[1] = -30. ;

//MESH SUBDIVISION_________________________________

//orthogonal to X
ElmsX = 30; 
//orthogonal to Y
ElmsY = 30; 

//orthogonal to Z
ElmsZ[1] = 30;  

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

// PHYSICAL VOLUMES_____________________________

Physical Volume("Layer_1",1) = {1};

// PHYSICAL SURFACES____________________________

Physical Surface("X_0",1) = {25};
Physical Surface("X_M",2) = {17};
Physical Surface("Y_0",3) = {13};
Physical Surface("Y_M",4) = {21};
Physical Surface("Z_0",5) = {1}; 

Mesh 3;
Save "cube.msh";