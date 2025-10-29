//Mesh.Format = 16;
Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

lc = 10;

// DISTANCES__________________________________________
// X
H = 100. ;
dx = 30. ;
// zona con edifici 30*13=390 ;

// Y
V = 100. ;
dy = 30. ;
// zona con edifici 30*22 = 660 ;

// COORDINATES________________________________________

// X coordinates
X[1] = 0. ;
X[2] = X[1] + H ;
For i In {3:15}
    X[i] = X[i-1] + dx ;
EndFor
X[16] = X[15] + H ; 

// Y coordinates
Y[1] = 0. ;
Y[2] = Y[1] + V ;
For i In {3:24}
    Y[i] = Y[i-1] + dy ;
EndFor
Y[25] = Y[24] + V ;

// Z coordinates
Z[0] =  0.; 


Z[1] = 15.47;
Z[2] = 3.33;
Z[3] = 7.62;
Z[4] = 2.38;
Z[5] = 2.15;
Z[6] = 4.05;  //B
Z[7] = 3.57; //B
Z[8] = 10.47; //B
Z[9] = 1.91; //B


//MESH SUBDIVISION_________________________________

//orthogonal to X
ElmsX[1] = 6.;  //3
ElmsX[2] = 3.;  //2

//orthogonal to Y
ElmsY[1] = 6.;  //3
ElmsY[2] = 3.;  //2

//orthogonal to Z
ElmsZ[1] = 2;
ElmsZ[2] = 2;
ElmsZ[3] = 2;
ElmsZ[4] = 2;
ElmsZ[5] = 2; 
ElmsZ[6] = 2; 
ElmsZ[7] = 2;  
ElmsZ[8] = 2;  
ElmsZ[9] = 2;

// POINTS________________________________________

For j In {1:25}
    For i In {1:16}
    p=newp;    Point(p) = { X[i], Y[j], Z[0], lc } ;
    EndFor
EndFor

// LINES__________________________________________

// Horizontal
For j In {0:24}
    For i In {1:15}
    l=newl;     Line(l) = { 16*j+i, 16*j+i+1 } ;
    EndFor
EndFor

// Vertical
For j In {0:23}
    For i In {1:16}
    l=newl;     Line(l) = { 16*j+i, 16*j+i+16 } ;
    EndFor
EndFor

// TRANSFINITE LINES_____________________________

// othogonal to X axis
For i In {0:24}
Transfinite Line{15*i+1} = ElmsX[1] + 1 ;
EndFor
For i In {1:25}
Transfinite Line{15*i} = ElmsX[1] + 1 ;
EndFor

For j In {2:14}
    For i In {0:24}
    Transfinite Line{j+15*i} = ElmsX[2] + 1;
    EndFor
EndFor

// othogonal to Y axis
For i In {376:391}
Transfinite Line{i} = ElmsY[1] + 1 ;
EndFor
For i In {744:759}
Transfinite Line{i} = ElmsY[1] + 1 ;
EndFor

For j In {392:728:16}
    For i In {0:15}
    Transfinite Line{i+j} = ElmsY[2] + 1 ;
    EndFor
EndFor


// CURVE LOOPS__________________________________

cl = 0 ;
For j In {0:23}
    For i In {1:15}
        cl = cl+1 ;
        Curve Loop(cl) = { 15*j+i, 16*j+i+376, -(15*j+i+15), -(16*j+i+375) } ;
        Plane Surface(cl) = cl ;
        Transfinite Surface(cl) ;
        Recombine Surface(cl) ;
    EndFor
EndFor


// EXTRUSION__________________________________

// Layer 1
Extrude{0,0,Z[1]} { Surface{1:360}; Layers{ElmsZ[1]}; Recombine;}

// Layer 2
Extrude{0,0,Z[2]} { Surface{781:8679:22}; Layers{ElmsZ[2]}; Recombine;}

// Layer 3
Extrude{0,0,Z[3]} { Surface{8701:16599:22}; Layers{ElmsZ[3]}; Recombine;}

// Layer 4
Extrude{0,0,Z[4]} { Surface{16621:24519:22}; Layers{ElmsZ[4]}; Recombine;}

// Layer 5
Extrude{0,0,Z[5]} { Surface{24541:32439:22}; Layers{ElmsZ[5]}; Recombine;}

// Layer 6
Extrude{0,0,Z[6]} { Surface{32461:40359:22}; Layers{ElmsZ[6]}; Recombine;}

// Layer 7
Extrude{0,0,Z[7]} { Surface{40381:48279:22}; Layers{ElmsZ[7]}; Recombine;}

// Layer 8
Extrude{0,0,Z[8]} { Surface{48301:56199:22}; Layers{ElmsZ[8]}; Recombine;}

// Layer 9
Extrude{0,0,Z[9]} { Surface{56221:64119:22}; Layers{ElmsZ[9]}; Recombine;}


// PHYSICAL VOLUMES___________________________

Physical Volume("Layer_1",1) = {1:360};
Physical Volume("Layer_2",2) = {361:720};
Physical Volume("Layer_3",3) = {721:1080};
Physical Volume("Layer_4",4) = {1081:1440};
Physical Volume("Layer_5",5) = {1441:1800};
Physical Volume("Layer_6",6) = {1801:2160}; 
Physical Volume("Layer_7",7) = {2161:2520}; 
Physical Volume("Layer_8",8) = {2521:2880}; 
Physical Volume("Layer_9",9) = {2881:3240}; 

// PHYSICAL SURFACES___________________________

Physical Surface("X_0",1) = {780:71730:330};
Physical Surface("X_M",2) = {1080:72030:330};
// Physical Surface("Y_0",3) = {768:1076:22, 8688:8996:22, 16608:16916:22, 24528:24836:22, 32448:32756:22 };
// Physical Surface("Y_M",4) = {8366:8674:22, 16286:16594:22, 24206:24514:22, 32126:32434:22, 40046:40354:22 };
Physical Surface("Z_0",5) = {1:360};
//Physical Surface("Z_M",6) = {32461:40359:22};

Mesh 3;
Save "mesh_noB.msh";