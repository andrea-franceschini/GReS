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
Extrude{0,0,Z[6]} {
Surface{ 
    32461:32791:22, 32835,32879,32923,32967,33011,33055, 33099:33473:22, 33517,33561,33605,33649,33693,
    33737:34111:22, 34155,34199,34243,34287,34331,34375, 34419:34793:22, 34837,34881,34925,34969,35013,
    35057:35431:22, 35475,35519,35563,35607,35651,35695, 35739:36113:22, 36157,36201,36245,36289,36333,
    36377:36751:22, 36795,36839,36883,36927,36971,37015, 37059:37433:22, 37477,37521,37565,37609,37653,
    37697:38071:22, 38115,38159,38203,38247,38291,38335, 38379:38753:22, 38797,38841,38885,38929,38973,
    39017:39391:22, 39435,39479,39523,39567,39611,39655, 39699:40359:22
}; Layers{ElmsZ[6]}; 
Recombine;} //B

// Layer 7
Extrude{0,0,Z[7]} { Surface{40381:46695:22}; Layers{ElmsZ[7]}; Recombine;} //B

// Layer 8
Extrude{0,0,Z[8]} { Surface{46717:53031:22}; Layers{ElmsZ[8]}; Recombine;} //B

// Layer 9
Extrude{0,0,Z[9]} { Surface{53053:59367:22}; Layers{ElmsZ[9]}; Recombine;} //B


// PHYSICAL VOLUMES___________________________

Physical Volume("Layer_1",1) = {1:360};
Physical Volume("Layer_2",2) = {361:720};
Physical Volume("Layer_3",3) = {721:1080};
Physical Volume("Layer_4",4) = {1081:1440};
Physical Volume("Layer_5",5) = {1441:1800};
Physical Volume("Layer_6",6) = {1801:2088}; //B
Physical Volume("Layer_7",7) = {2089:2376}; //B
Physical Volume("Layer_8",8) = {2377:2664}; //B
Physical Volume("Layer_9",9) = {2665:2952}; //B

// PHYSICAL SURFACES___________________________

Physical Surface("X_0",1) = {780:40050:330, 
    40380,46716,53052,59388, 46386,52722,59058,65394,
    40710:45880:1034, 40886:46056:1034, 41216:45352:1034, 41414:45550:1034,
    47046:52216:1034, 47222:52392:1034, 47552:51688:1034, 47750:51886:1034,
    53382:58552:1034, 53558:58728:1034, 53888:58024:1034, 54086:58222:1034,
    59718:64888:1034, 59894:65064:1034, 60224:64360:1034, 60422:64558:1034 
};
Physical Surface("X_M",2) = {1080:40350:330,
    40680,47016,53352,59688, 46686,53022,59358,65694,
    40856:46026:1034, 41186:46356:1034, 41384:45520:1034, 41714:45850:1034,
    47192:52362:1034, 47522:52692:1034, 47720:51856:1034, 48050:52186:1034,
    53528:58698:1034, 53858:59028:1034, 54056:58192:1034, 54386:58522:1034,
    59864:65034:1034, 60194:65364:1034, 60392:64528:1034, 60722:64858:1034
};
// Physical Surface("Y_0",3) = {};
// Physical Surface("Y_M",4) = {};
Physical Surface("Z_0",5) = {1:360};
//Physical Surface("Z_M",6) = {};

Mesh 3;
Save "Mesh_B2.msh";

