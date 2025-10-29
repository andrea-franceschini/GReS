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
Extrude {0, 0, Z[6]} {
  Surface{ 
    32461, 32483, 32505, 32527, 32549, 32571, 32593, 32615, 32637, 32659, 32681, 32703, 32725, 32747, 32769,
    32791, 32857, 32879, 32901, 32967, 32989, 33011, 33033, 33077, 33099, 33121, 33143, 33209, 33231, 33253, 33275, 33297, 33341, 33363, 33429,
    33451, 33473, 33495, 33517, 33539, 33605, 33627, 33671, 33737, 33759, 33781, 33803, 33869, 33891, 33957, 33979, 34001, 34045, 34067, 34089,
    34111, 34177, 34199, 34221, 34243, 34331, 34353, 34375, 34397, 34419, 34441, 34485, 34507, 34573, 34595, 34617, 34639, 34661, 34683, 34749,
    34771, 34793, 34815, 34881, 34903, 34969, 35013, 35035, 35057, 35079, 35101, 35145, 35167, 35189, 35211, 35233, 35277, 35299, 35365, 35409,
    35431, 35519, 35563, 35585, 35607, 35629, 35651, 35673, 35695, 35739, 35761, 35805, 35827, 35849, 35871, 35893, 35959, 36025, 36047, 36069,
    36091, 36113, 36135, 36201, 36223, 36245, 36267, 36289, 36311, 36399, 36421, 36465, 36487, 36553, 36575, 36641, 36663, 36685, 36707, 36729,
    36751, 36817, 36839, 36861, 36883, 36905, 36927, 36971, 36993, 37059, 37081, 37103, 37125, 37147, 37169, 37235, 37301, 37323, 37345, 37389,
    37411, 37433, 37499, 37543, 37565, 37587, 37609, 37631, 37653, 37719, 37741, 37763, 37807, 37829, 37873, 37917, 37983, 38005, 38027, 38049,
    38071, 38093, 38159, 38181, 38203, 38225, 38247, 38291, 38313, 38379, 38401, 38423, 38445, 38467, 38489, 38577, 38599, 38621, 38687, 38709,
    38731, 38797, 38819, 38863, 38907, 38951, 38973, 38995, 39017, 39039, 39061, 39149, 39171, 39193, 39215, 39237, 39259, 39281, 39325, 39369,
    39391, 39413, 39435, 39457, 39479, 39501, 39567, 39589, 39611, 39699, 39721, 39743, 39809, 39875, 39897, 39941, 39963, 39985, 40007, 40029,
    40051, 40073, 40095, 40117, 40139, 40161, 40183, 40205, 40227, 40249, 40271, 40293, 40315, 40337, 40359
  };
  Layers{ElmsZ[6]};
  Recombine;
}

// Layer 7
Extrude{0,0,Z[7]} { Surface{40381:45859:22}; Layers{ElmsZ[7]}; Recombine;}

// Layer 8
Extrude{0,0,Z[8]} { Surface{45881:51359:22}; Layers{ElmsZ[8]}; Recombine;}

// Layer 9
Extrude{0,0,Z[9]} { Surface{51381:56859:22}; Layers{ElmsZ[9]}; Recombine;}


// PHYSICAL VOLUMES___________________________

Physical Volume("Layer_1",1) = {1:360};
Physical Volume("Layer_2",2) = {361:720};
Physical Volume("Layer_3",3) = {721:1080};
Physical Volume("Layer_4",4) = {1081:1440};
Physical Volume("Layer_5",5) = {1441:1800};
Physical Volume("Layer_6",6) = {1801:2050}; //B
Physical Volume("Layer_7",7) = {2051:2300}; //B
Physical Volume("Layer_8",8) = {2301:2550}; //B
Physical Volume("Layer_9",9) = {2551:2800}; //B

// PHYSICAL SURFACES___________________________

Physical Surface("X_0",1) = {
    780:40050:330, 40380,45880,51380,56880, 40710:45330:220, 46210:50830:220, 51710:56330:220, 57210:61830:220, 45550,51050,56550,62050      
};       

Physical Surface("X_M",2) = {
    1080:40350:330,  40680,46180,51680,57180, 40900:45520:220, 46400:51020:220, 51900:56520:220, 57400:62020:220, 45850,51350,56850,62350
};

// Physical Surface("Y_0",3) = {};
// Physical Surface("Y_M",4) = {};
Physical Surface("Z_0",5) = {1:360};
// Physical Surface("Z_M",6) = {};

Mesh 3;
Save "Mesh_B3.msh";




