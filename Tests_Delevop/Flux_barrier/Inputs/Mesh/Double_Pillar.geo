Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

lc = 1e-1;

// Distances
// X
dx = 30. ;
Ha = 15. ;
dxb = 5. ; 
Hb = 15. ;
// Y
dy = 30. ;
Va = 15. ;
dyb = 5. ;
Vb = 15. ;

// Elevation for extrusion
Z[1] = 6. ;       Z[4] = 5. ;    Z[7] = 15. ;
Z[2] = 6. ;       Z[5] = 5. ;    
Z[3] = 6. ;       Z[6] = 5. ;  

// Values for mesh subdivision
ElmsX[1] = 10. ;
ElmsX[2] = 5. ;
ElmsX[3] = 5. ;

ElmsY[1] = 10.;
ElmsY[2] = 5.;
ElmsY[3] = 10.;

ElmsZ[1] = 5.;
ElmsZ[2] = 5.;
ElmsZ[3] = 5.;
//___________________________________________________________
// Outer points dimensions
X[1] = 0. ; 
X[2] = X[1] + dx ;      
X[3] = X[2] + Ha ;
X[4] = X[3] + dxb ;
X[5] = X[4] + Hb ;
X[6] = X[5] + dx ;

Y[1] = 0. ;
Y[2] = Y[1] + dy ;
Y[3] = Y[2] + Va ;
Y[4] = Y[3] + dyb ;
Y[5] = Y[4] + Vb ;
Y[6] = Y[5] + dy ;

Z[0] = 0 ;  

// Outer Points
Point(1) = { X[1], Y[1], Z[0], lc };
Point(2) = { X[2], Y[1], Z[0], lc };
Point(3) = { X[3], Y[1], Z[0], lc };
Point(4) = { X[4], Y[1], Z[0], lc };
Point(5) = { X[5], Y[1], Z[0], lc };
Point(6) = { X[6], Y[1], Z[0], lc };

Point(7) = { X[6], Y[2], Z[0], lc };
Point(8) = { X[6], Y[3], Z[0], lc };
Point(9) = { X[6], Y[4], Z[0], lc };
Point(10) = { X[6], Y[5], Z[0], lc };
Point(11) = { X[6], Y[6], Z[0], lc };

Point(12) = { X[5], Y[6], Z[0], lc };
Point(13) = { X[4], Y[6], Z[0], lc };
Point(14) = { X[3], Y[6], Z[0], lc };
Point(15) = { X[2], Y[6], Z[0], lc };
Point(16) = { X[1], Y[6], Z[0], lc };

Point(17) = { X[1], Y[5], Z[0], lc };
Point(18) = { X[1], Y[4], Z[0], lc };
Point(19) = { X[1], Y[3], Z[0], lc };
Point(20) = { X[1], Y[2], Z[0], lc };

// Inside Points
// building
Point(21) = { X[2], Y[2], Z[0], lc };      Point(25) = { X[4], Y[4], Z[0], lc };
Point(22) = { X[3], Y[2], Z[0], lc };      Point(26) = { X[5], Y[4], Z[0], lc };
Point(23) = { X[3], Y[3], Z[0], lc };      Point(27) = { X[5], Y[5], Z[0], lc };
Point(24) = { X[2], Y[3], Z[0], lc };      Point(28) = { X[4], Y[5], Z[0], lc };
// line intersection
Point(29) = { X[4], Y[2], Z[0], lc };     Point(33) = { X[2], Y[4], Z[0], lc };
Point(30) = { X[5], Y[2], Z[0], lc };     Point(34) = { X[3], Y[4], Z[0], lc };
Point(31) = { X[5], Y[3], Z[0], lc };     Point(35) = { X[3], Y[5], Z[0], lc };
Point(32) = { X[4], Y[3], Z[0], lc };     Point(36) = { X[2], Y[5], Z[0], lc };

// Lines
// Contour lines
For s In {1:19}
    Line(s) = {s, s+1} ;
EndFor
Line(20) = {20,1} ;
//------------------
For s In {21:23}
    Line(s) = {s, s+1} ;
EndFor
Line(24) = {24,21} ;
//-------------------
For s In {25:27}
    Line(s) = {s, s+1} ;
EndFor
Line(28) = {28,25} ;

// Inside horizontal lines
Line(29) = {20,21} ;        Line(37) = {18,33} ;
Line(30) = {22,29} ;        Line(38) = {33,34} ;
Line(31) = {29,30} ;        Line(39) = {34,25} ;
Line(32) = {30,7} ;         Line(40) = {26,9} ;
Line(33) = {19,24} ;        Line(41) = {17,36} ;
Line(34) = {23,32} ;        Line(42) = {36,35} ;
Line(35) = {32,31} ;        Line(43) = {35,28} ;
Line(36) = {31,8} ;         Line(44) = {27,10} ;

// Inside vertical lines
Line(45) = {2,21} ;         Line(53) = {4,29} ;
Line(46) = {24,33} ;        Line(54) = {29,32} ;
Line(47) = {33,36} ;        Line(55) = {32,25} ;
Line(48) = {36,15} ;        Line(56) = {28,13} ;
Line(49) = {3,22} ;         Line(57) = {5,30} ;
Line(50) = {23,34} ;        Line(58) = {30,31} ;
Line(51) = {34,35} ;        Line(59) = {31,26} ;
Line(52) = {35,14} ;        Line(60) = {27,12} ;

// Mesh subdivision
// othogonal to X axis
Transfinite Line{20,45,49,53,57,6} = ElmsX[1] + 1 ;
Transfinite Line{19,24,22,54,58,7} = ElmsX[2] + 1 ;
Transfinite Line{18,46,50,55,59,8} = ElmsX[2] + 1 ;
Transfinite Line{17,47,51,28,26,9} = ElmsX[2] + 1 ;
Transfinite Line{16,48,52,56,60,10} = ElmsX[1] + 1 ;
// othogonal to Y axis
Transfinite Line{1,29,33,37,41,15} = ElmsY[1] + 1 ;
Transfinite Line{2,21,23,38,42,14} = ElmsY[1] + 1 ;
Transfinite Line{3,30,34,39,43,13} = ElmsY[2] + 1 ;
Transfinite Line{4,31,35,25,27,12} = ElmsY[1] + 1 ;
Transfinite Line{5,32,36,40,44,11} = ElmsY[1] + 1 ;

// Delimitation of surfaces
Curve Loop(1) = {1,45,-29,20};     Curve Loop(6) = {29,-24,-33,19};
Curve Loop(2) = {2,49,-21,-45};    Curve Loop(7) = {21,22,23,24};
Curve Loop(3) = {3,53,-30,-49};    Curve Loop(8) = {30,54,-34,-22};
Curve Loop(4) = {4,57,-31,-53};    Curve Loop(9) = {31,58,-35,-54};
Curve Loop(5) = {5,6,-32,-57};     Curve Loop(10) = {32,7,-36,-58};

Curve Loop(11) = {33,46,-37,18};     Curve Loop(16) = {37,47,-41,17};
Curve Loop(12) = {-23,50,-38,-46};   Curve Loop(17) = {38,51,-42,-47};
Curve Loop(13) = {34,55,-39,-50};    Curve Loop(18) = {39,-28,-43,-51};
Curve Loop(14) = {35,59,-25,-55};    Curve Loop(19) = {25,26,27,28};
Curve Loop(15) = {36,8,-40,-59};     Curve Loop(20) = {40,9,-44,-26};

Curve Loop(21) = {41,48,15,16};
Curve Loop(22) = {42,52,14,-48};
Curve Loop(23) = {43,56,13,-52};
Curve Loop(24) = {-27,60,12,-56};
Curve Loop(25) = {44,10,11,-60};

// Creating Surfaces
For s In {1:25}
    Plane Surface(s) = {s};
    Transfinite Surface{s};
    Recombine Surface{s};
EndFor

// Extruding surfaces
// soil below the column (layer 1-2-3)
Extrude{0,0,Z(1)} {Surface{1:25}; Layers{ElmsZ[1]}; Recombine;}
Extrude{0,0,Z(2)} {Surface{82,104,126,148,170,192,214,236,258,280,302,324,346,368,390,412,434,456,478,500,522,544,566,588,610}; Layers{ElmsZ[1]}; Recombine;}
Extrude{0,0,Z(3)} {Surface{632,654,676,698,720,742,764,786,808,830,852,874,896,918,940,962,984,1006,1028,1050,1072,1094,1116,1138,1160}; Layers{ElmsZ[1]}; Recombine;}

// soil around column + columns (layer 4)
Extrude{0,0,Z(4)} {Surface{1182,1204,1226,1248,1270,1292,1314,1336,1358,1380,1402,1424,1446,1468,1490,1512,1534,1556,1578,1600,1622,1644,1666,1688,1710}; Layers{ElmsZ[2]}; Recombine;}

// soil around column + column layer 5
Extrude{0,0,Z(5)} {Surface{1732,1754,1776,1798,1820,1842,1864,1886,1908,1930,1952,1974,1996,2018,2040,2062,2084,2106,2128,2150,2172,2194,2216,2238,2260}; Layers{ElmsZ[2]}; Recombine;}

// soil around column + column layer 6
Extrude{0,0,Z(6)} {Surface{2282,2304,2326,2348,2370,2392,2414,2436,2458,2480,2502,2524,2546,2568,2590,2612,2634,2656,2678,2700,2722,2744,2766,2788,2810}; Layers{ElmsZ[2]}; Recombine;}

// column above the ground
Extrude{0,0,Z(7)} {Surface{2964,3228}; Layers{ElmsZ[3]}; Recombine;}

// PHYSICAL VOLUMES
// physical volume of soil below and around the column
Physical Volume("Soil_layer_1", 1) = {1:25};
Physical Volume("Soil_layer_2", 2) = {26:50};
Physical Volume("Soil_layer_3", 3) = {51:75};
Physical Volume("Soil_layer_4", 4) = {76:81,83:93,95:100};
Physical Volume("Soil_layer_5", 5) = {101:106,108:118,120:125};
Physical Volume("Soil_layer_6", 6) = {126:131,133:143,145:150};

// physical volume of the columns
Physical Volume("building_A",7) = {82,107,132,151};
Physical Volume("building_B",8) = {94,119,144,152};

// Physical surfaces
Physical Surface("X0",1) = {81,191,301,411,521,631,741,851,961,1071,1181,1291,1401,1511,1621,1731,1841,1951,2061,2171,2281,2391,2501,2611,2721,2831,2941,3051,3161,3271};
Physical Surface("Xm",2) = {161,271,381,491,601,711,821,931,1041,1151,1261,1371,1481,1591,1701,1811,1921,2031,2141,2251,2361,2471,2581,2691,2801,2911,3021,3131,3241,3351};
Physical Surface("Y0",3) = {69,91,113,135,157,619,641,663,685,707,1169,1191,1213,1235,1257,1719,1741,1763,1785,1807,2269,2291,2313,2335,2357,2819,2841,2863,2885,2907};
Physical Surface("Ym",4) = {517,539,561,583,605,1067,1089,1111,1133,1155,1617,1639,1661,1683,1705,2167,2189,2211,2233,2255,2717,2739,2761,2783,2805,3267,3289,3311,3333,3355};
Physical Surface("Z0",5) = {1:25};
Physical Surface("Zm",6) = {2832,2854,2876,2898,2920,2942,2986,3008,3030,3052,3074,3096,3118,3140,3162,3184,3206,3250,3272,3294,3316,3338,3360};

// Surfaces of building - A
Physical Surface("X0_A",7) = {1833,2383,2933,3381};
Physical Surface("Xm_A",8) = {1855,2405,2955,3373};
Physical Surface("Y0_A",9) = {1749,2299,2849,3369};
Physical Surface("Ym_A",10) = {1859,2409,2959,3377};
Physical Surface("Z0_A",11) = {1314};
Physical Surface("Zm_A",12) = {3382};

// Surfaces of building - B
Physical Surface("X0_B",13) = {2097,2647,3197,3403};
Physical Surface("Xm_B",14) = {2119,2669,3219,3395};
Physical Surface("Y0_B",15) = {2013,2563,3113,3391};
Physical Surface("Ym_B",16) = {2123,2673,3223,3399};
Physical Surface("Z0_B",17) = {1578};
Physical Surface("Zm_B",18) = {3404};

Mesh 3 ;
Save "Double_Pillar.msh" ;