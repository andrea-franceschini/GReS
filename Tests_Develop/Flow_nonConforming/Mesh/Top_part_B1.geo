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
Extrude {0, 0, Z[4]} {
  Surface{
    1:16,
    20:31,36,37,38,41:51,53:56,59:62,65,66,69,70,73:77,79:82,85:92,
    95,97,98,100,101,103,105:112,115:118,120,121,124,127:130,133:140,143,145:152,155:158,
    161,164:168,170,171,174:183,186,187,190,192,193,195,196,198,199,200,201,204:211,213,214,215,218,219,221,
    224:228,230,232,233,234,236:241,243,244,245,247,249,251,254:258,260,261,262,264,265,266,268:271,273,274,275,277,
    278,279,281,285:288,290,291,292,294:298,300,301,303,305,307,309,310,313:316,318,319,320,322:328,330,331,333,335,
    339,340,341,343:360
  };
  Layers{ElmsZ[4]};
  Recombine;
}

// Layer 2
Extrude{0,0,Z[5]} { Surface{781:6501:22}; Layers{ElmsZ[5]}; Recombine;}

// Layer 3
Extrude{0,0,Z[6]} { Surface{6523:12243:22}; Layers{ElmsZ[6]}; Recombine;}

// Layer 4
Extrude{0,0,Z[7]} { Surface{12265:17985:22}; Layers{ElmsZ[7]}; Recombine;}

// Layer 5
Extrude{0,0,Z[8]} { Surface{18007:23727:22}; Layers{ElmsZ[8]}; Recombine;}

// Layer 6
Extrude{0,0,Z[9]} { Surface{23749:29469:22}; Layers{ElmsZ[9]}; Recombine;}

// PHYSICAL VOLUMES___________________________

Physical Volume("Layer_4",1) = {1:261};
Physical Volume("Layer_5",2) = {262:522};
Physical Volume("Layer_6",3) = {523:783};
Physical Volume("Layer_7",4) = {784:1044};
Physical Volume("Layer_8",5) = {1045:1305};
Physical Volume("Layer_9",6) = {1306:1566};

// PHYSICAL SURFACES___________________________

Physical Surface("X_0",1) = {780,6522,12264,18006,23748,29490, 6192,11934,17676,23418,29160,34902, 1110:5730:462, 6852:11736:462 , 12594:17478:462 ,18336:23220:462 ,
    24078:28962:462, 29820:34704:462 , 1374:5994:462 , 7116:11736:462, 12858:17478:462 , 18600:23220:462, 24342:28962:462, 30084:34704:462  };


Physical Surface("X_M",2) = {1080,6822,12564,18306,24048,29790,6492,12234,17976,23718,29460,35202, 1344:5964:462,7086:11706:462,12828:17448:462,18570:23190:462,
    24312:28932:462,30054:34674:462, 1542:6162:462, 7284:11904:462, 13026:17646:462, 18768:23388:462, 24510:29130:462, 30252:34872:462 };
// Physical Surface("Y_0",3) = {768:1076:22, 8688:8996:22, 16608:16916:22, 24528:24836:22, 32448:32756:22, 40368:40676:22 };
// Physical Surface("Y_M",4) = {8366:8674:22, 16286:16594:22, 24206:24514:22, 32126:32434:22, 40046:40354:22, 47966:48274:22};
Physical Surface("Z_0",5) = {1:360};    // BOTTOM INTERFACE
// Physical Surface("Z_M",6) = {40381:48279:22};

Mesh 3;
Save "Top_part_B1.msh";
//Save "Top_part.vtk";