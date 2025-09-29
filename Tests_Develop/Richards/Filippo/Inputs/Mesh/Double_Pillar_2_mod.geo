Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

lc = 1e-1;

// Distances
// X
dx = 120. ;
Ha = 20. ;
Hb = 20. ;
I = 20. ; 

// Y
dy = 95. ;
Va = 40. ;
dyb = 30. ;
Vb = 40. ;

nregs=6;  //6, 9

// Outer points dimensions
X[1] = 0. ; 
X[2] = X[1] + dx ;      
X[3] = X[2] + Ha ;
X[4] = X[3] + I ;
X[5] = X[4] + Hb ; 
X[6] = X[5] + dx ;

// X[6] = X[5] + Hb ;
// X[7] = X[6] + Hb ; 
// X[8] = X[7] + Hb ; 
// X[9] = X[8] + dx ;  

Y[1] = 0. ;
Y[2] = Y[1] + dy ;
Y[3] = Y[2] + Va ;
Y[4] = Y[3] + dyb ;
Y[5] = Y[4] + Vb ;
Y[6] = Y[5] + dy ;

// Y[6] = Y[5] + Vb ;
// Y[7] = Y[6] + Vb ;
// Y[8] = Y[7] + Vb ;
// Y[9] = Y[8] + dy ;

Z[0] = 0 ;  

// Elevation for extrusion
Z[1] = 10. ;       Z[4] = 10. ;
Z[2] = 10. ;       Z[5] = 10. ;    
Z[3] = 10. ;       Z[6] = 10. ;  

// Values for mesh subdivision
ElmsX[1] = 8;
ElmsX[2] = 4;
ElmsX[3] = 4;

ElmsY[1] = 8;
ElmsY[2] = 4;
ElmsY[3] = 4;

ElmsZ[1] = 4;
ElmsZ[2] = 6;
ElmsZ[3] = 8;
ElmsZ[4] = 20;

// Outer Points
For j In {1:nregs}
    For i In {1:nregs}
        pt=newp;
        Point(pt) = { X[i], Y[j], Z[0], lc };
    EndFor
EndFor

// Lines
// in X 
For j In {0:nregs-1}
    For i In {1:nregs-1}
        ln=newl;
        Line(ln) = {nregs*j+i, nregs*j+i+1} ;
    EndFor
EndFor

// in Y 
For j In {0:nregs-2}
    For i In {1:nregs}
        ln=newl;
        Line(ln) = {nregs*j+i, nregs*j+i+nregs} ;
    EndFor
EndFor

// Mesh subdivision
// othogonal to X axis
Transfinite Line{1, 6,11,16,21,26} = ElmsX[1]+1;
Transfinite Line{2, 7,12,17,22,27} = ElmsX[2]+1;
Transfinite Line{3, 8,13,18,23,28} = ElmsX[3]+1;
Transfinite Line{4, 9,14,19,24,29} = ElmsX[2]+1;
Transfinite Line{5,10,15,20,25,30} = ElmsX[1]+1;

// othogonal to Y axis
Transfinite Line{31:36} = ElmsY[1]+1;
Transfinite Line{37:42} = ElmsY[2]+1;
Transfinite Line{43:48} = ElmsY[3]+1;
Transfinite Line{49:54} = ElmsY[2]+1;
Transfinite Line{55:60} = ElmsY[1]+1;

// Delimitation of surfaces
count=0;
For j In {0:4}
    For i In {1:5}
        count=count+1;
        Curve Loop(count) = {(5*j+i),(6*j+i+1)+30,-(5*j+i+5),-(6*j+i)-30};
        Plane Surface(count) = {count};
        Transfinite Surface{count};
        Recombine Surface{count};
    EndFor
EndFor

// Extruding surfaces
// soil below the column (layer 1-2-3)
Extrude{0,0,Z(1)} {Surface{1:25}; Layers{ElmsZ[1]}; Recombine;}
Extrude{0,0,Z(2)} {Surface{82,104,126,148,170,192,214,236,258,280,302,324,346,368,390,412,434,456,478,500,522,544,566,588,610}; Layers{ElmsZ[1]}; Recombine;}
Extrude{0,0,Z(3)} {Surface{632,654,676,698,720,742,764,786,808,830,852,874,896,918,940,962,984,1006,1028,1050,1072,1094,1116,1138,1160}; Layers{ElmsZ[1]}; Recombine;}
// physical volume of soil below and around the column
Physical Volume("Soil_layer_1", 1) = {1:25};
Physical Volume("Soil_layer_2", 2) = {26:50};
Physical Volume("Soil_layer_3", 3) = {51:75};

// soil around column (layer 4)
Extrude{0,0,Z(4)} {Surface{1182,1204,1226,1248,1270,1292,1336,1358,1380,1402,1424,1446,1468,1490,1512,1534,1556,1600,1622,1644,1666,1688,1710}; Layers{ElmsZ[2]}; Recombine;}
Physical Volume("Soil_layer_4", 4) = {76:98};

// soil around column (layer 5)
Extrude{0,0,Z(5)} {Surface{1732,1754,1776,1798,1820,1842,1864,1886,1908,1930,1952,1974,1996,2018,2040,2062,2084,2106,2128,2150,2172,2194,2216}; Layers{ElmsZ[3]}; Recombine;}
Physical Volume("Soil_layer_5", 5) = {99:121};

// soil around column + column layer 6
Extrude{0,0,Z(6)} {Surface{2238,2260,2282,2304,2326,2348,2370,2392,2414,2436,2458,2480,2502,2524,2546,2568,2590,2612,2634,2656,2678,2700,2722}; Layers{ElmsZ[4]}; Recombine;}
Physical Volume("Soil_layer_6", 6) = {122:144};

// Physical surfaces
Physical Surface("X0",1) = {81,191,301,411,521,631,741,851,961,1071,1181,1291,1401,1511,1621,1731,1841,1929,2039,2127,2237,2347,2435,2545,2633,2743,2853,2941,3051,3139};
Physical Surface("Xm",2) = {161,271,381,491,601,711,821,931,1041,1151,1261,1371,1481,1591,1701,1811,1899,2009,2097,2207,2317,2405,2515,2603,2713,2823,2911,3021,3109,3219};
Physical Surface("Y0",3) = {69,91,113,135,157,619,641,663,685,707,1169,1191,1213,1235,1257,1719,1741,1763,1785,1807,2225,2247,2269,2291,2313,2731,2753,2775,2797,2819};
Physical Surface("Ym",4) = {517,539,561,583,605,1067,1089,1111,1133,1155,1617,1639,1661,1683,1705,2123,2145,2167,2189,2211,2629,2651,2673,2695,2717,3135,3157,3179,3201,3223};
Physical Surface("Z0",5) = {1:25};
Physical Surface("Zm",6) = {2744,2766,2788,2810,2832,2854,2876,2898,2920,2942,2964,2986,3008,3030,3052,3074,3096,3118,3140,3162,3184,3206,3228};

/*
// Surfaces of building - A
Physical Surface("X0_A",7) = {1833,2383,2933,3381};
Physical Surface("Xm_A",8) = {1877,2427,2977,3395};
Physical Surface("Y0_A",9) = {1749,1771,2299,2321,2849,2871,3369,3391};
Physical Surface("Ym_A",10) = {1881,1859,2431,2409,2981,2959,3377,3399};
Physical Surface("Z0_A",11) = {1314,1336};
Physical Surface("Zm_A",12) = {3382,3404};

// Surfaces of building - B
Physical Surface("X0_B",13) = {2075,2625,3175,3425};
Physical Surface("Xm_B",14) = {2119,2669,3219,3439};
Physical Surface("Y0_B",15) = {1991,2013,2541,2563,3091,3113,3413,3435};
Physical Surface("Ym_B",16) = {2123,2101,2673,2651,3223,3201,3443,3421};
Physical Surface("Z0_B",17) = {1556,1578};
Physical Surface("Zm_B",18) = {3426,3448};
*/

Mesh 3 ;
Save "Double_Pillar_2.msh" ;