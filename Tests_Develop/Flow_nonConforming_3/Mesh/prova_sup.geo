// Mesh.Format = 16; //  output format 
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
Z[1] = 50 ; 

//Extrusion hight

E = 40. ;

//MESH SUBDIVISION______________

//orthogonal to X
ElmsX[1] = 20. ; //5  //40  //10

//orthogonal to Y
ElmsY[1] = 20. ; //4  //30  //15

//orthogonal to Z
ElmsZ[1] = 10. ;  //8  //15  //15


// POINTS______________________
 
p1 = newp;  Point(p1) = { X[1], Y[1], Z[1], lc } ;
p2 = newp;  Point(p2) = { X[2], Y[2], Z[1], lc } ;
p3 = newp;  Point(p3) = { X[3], Y[3], Z[1], lc } ;
p4 = newp;  Point(p4) = { X[4], Y[4], Z[1], lc } ;

// LINES________________________

l1 = newl;  Line(l1) = {p1,p2} ;
l2 = newl;  Line(l2) = {p2,p3} ;
l3 = newl;  Line(l3) = {p3,p4} ;
l4 = newl;  Line(l4) = {p4,p1} ;

// TRANSFINITE LINES_____________

// orthogonal to X
Transfinite Line{l1,l3} = ElmsX[1]+1 ;

// orthogonal to Y
Transfinite Line{l2,l4} = ElmsY[1]+1 ;

// SURFACE_______________________

cl = newcl;
Curve Loop(cl) = {l1,l2,l3,l4} ;

s = news;
Plane Surface(s) = cl ;

Transfinite Surface{s} ;
Recombine Surface{s} ;

// EXTRUDE________________________

// LAYER top 
toplayer[] = Extrude{0,0,E} {Surface{s}; Layers{ElmsZ[1]}; Recombine; } ;

//________________________________


Physical Volume("toplayer",1) = toplayer[1] ;

Physical Surface("Z_M",1) = toplayer[0] ;
Physical Surface("X_0",2) = toplayer[5] ;
Physical Surface("X_M",3) = toplayer[3] ;
Physical Surface("Y_0",4) = toplayer[2] ;
Physical Surface("Y_M",5) = toplayer[4] ;
Physical Surface("topinterface",6) = s ;



Mesh 3 ;
// Save "prova_inf.vtk" ;
Save "prova_sup.msh" ;