// Mesh.Format = 16; // vtk output format
Mesh.Format = 1;
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

lc=1 ;

// X
X[1] = 0. ;
X[2] = 200. ;
X[3] = 200. ;
X[4] = 0. ;

// Y
Y[1] = 0. ;
Y[2] = 0. ;
Y[3] = 200. ;
Y[4] = 200. ;

// Z
// Z[1] = 20 ; 

// extrusion hight

E = 20. ;

//MESH SUBDIVISION______________

// orthogonal to X
ElmsX[1] = 5. ;

// orthogonal to Y
ElmsY[1] = 4. ;

//orthogonal to Z
ElmsZ[1] = 4. ;


// POINTS______________________

// SURFACE  
Point(1) = {X[1],Y[1],0,lc} ;
Point(2) = {X[2],Y[2],0,lc} ;
Point(3) = {X[3],Y[3],0,lc} ;
Point(4) = {X[4],Y[4],0,lc} ;

// LINES___________________

Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;

// TRANSFINITE LINES_____________

// orthogonal to X
Transfinite Line{1,3} = ElmsX[1]+1 ;

// orthogonal to Y
Transfinite Line{2,4} = ElmsY[1]+1 ;

// SURFACE_______________________

Curve Loop(1) = {1,2,3,4} ;
Plane Surface(1) = 1 ;

Transfinite Surface{1} ;
Recombine Surface{1} ;

// EXTRUDE________________________

// LAYER 1
Layer_1[] = Extrude{0,0,E} {Surface{1}; Layers{ElmsZ[1]}; Recombine;};

//________________________________

Physical Volume("Layer_1",1) = {1} ;

Physical Surface("bottom_interface",1) = Layer_1[0] ;
Physical Surface("X_0",2) = Layer_1[5] ;
Physical Surface("X_M",3) = Layer_1[3] ;
Physical Surface("Y_0",4) = Layer_1[2] ;
Physical Surface("Y_M",5) = Layer_1[4] ;
Physical Surface("Z_0",6) = 1 ;


Mesh 3 ;
// Save "prova_inf.vtk" ;
Save "prova_inf.msh" ;