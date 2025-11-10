Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
// Mesh.Format = 16; // vtk output format

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:
lc = 10;
fileName = "Top_part.msh";
// fileName = "Top_part.vtk";

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
Z[0] = 21.03 ; 

//Z[1] = 15.48;
//Z[2] = 3.33;  
//Z[3] = 2.22;   
Z[4] = 5.65; 
Z[5] = 2.14; 
Z[6] = 6.19; 
Z[7] = 3.57;    
Z[8] = 10.48;
Z[9] = 1.90;


//MESH SUBDIVISION_________________________________

//orthogonal to X
ElmsX[1] = 3.;  //8
ElmsX[2] = 2.;  //5

//orthogonal to Y
ElmsY[1] = 3.;  //8
ElmsY[2] = 2.;  //5

//orthogonal to Z
// ElmsZ[1] = 31;
// ElmsZ[2] = 7;
// ElmsZ[3] = 4;
ElmsZ[4] = 3;  // 11
ElmsZ[5] = 3;   // 4 
ElmsZ[6] = 3;  // 12
ElmsZ[7] = 3;   // 7
ElmsZ[8] = 5;  // 21
ElmsZ[9] = 3;   // 4

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
Extrude{0,0,Z[4]} { Surface{1:360}; Layers{ElmsZ[4]}; Recombine;}

// Layer 2
Extrude{0,0,Z[5]} { Surface{781:8679:22}; Layers{ElmsZ[5]}; Recombine;}

// Layer 3
Extrude{0,0,Z[6]} { Surface{8701:16599:22}; Layers{ElmsZ[6]}; Recombine;}

// Layer 4
Extrude{0,0,Z[7]} { Surface{16621:24519:22}; Layers{ElmsZ[7]}; Recombine;}

// Layer 5
Extrude{0,0,Z[8]} { Surface{24541:32439:22}; Layers{ElmsZ[8]}; Recombine;}

// Layer 6
Extrude{0,0,Z[9]} { Surface{32461:40359:22}; Layers{ElmsZ[9]}; Recombine;}

// PHYSICAL VOLUMES___________________________

Physical Volume("Layer_4",1) = {1:360};
Physical Volume("Layer_5",2) = {361:720};
Physical Volume("Layer_6",3) = {721:1080};
Physical Volume("Layer_7",4) = {1081:1440};
Physical Volume("Layer_8",5) = {1441:1800};
Physical Volume("Layer_9",6) = {1801:2160};

// PHYSICAL SURFACES___________________________

Physical Surface("X_0",1) = {780:47970:330};
Physical Surface("X_M",2) = {1080:48270:330};
Physical Surface("Y_0",3) = {768:1076:22, 8688:8996:22, 16608:16916:22, 24528:24836:22, 32448:32756:22, 40368:40676:22 };
Physical Surface("Y_M",4) = {8366:8674:22, 16286:16594:22, 24206:24514:22, 32126:32434:22, 40046:40354:22, 47966:48274:22};
Physical Surface("Z_0",5) = {1:360};    // BOTTOM INTERFACE
Physical Surface("Z_M",6) = {40381:48279:22};

Mesh 3;
Save Sprintf(fileName);