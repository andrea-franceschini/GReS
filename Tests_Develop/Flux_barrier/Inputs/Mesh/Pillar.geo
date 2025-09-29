Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

lc = 1e-2;

// Domain dimension
B = 50. ;
b = 5. ;
H = 50. ;
h = 5. ;

// Plane XY
X[1] = 0.; 
X[2] = (B-b)/2 ;
X[3] = X[2] + b ;
X[4] = B ;

Y[1] = 0. ;
Y[2] = (H-h)/2 ;
Y[3] = Y[2] + h ;
Y[4] = H ;

// Dimension Z
Z[1] = 5. ;
Z[2] = 5. ;
Z[3] = 13. ;

// Values for mesh subdivision
ElmsX[1] = 10 ; // 10
ElmsX[2] = 5 ;  // 5
ElmsX[3] = 10 ; // 10

ElmsY[1] = 10 ;  // 10
ElmsY[2] = 5 ;   // 5
ElmsY[3] = 10 ;  // 10

ElmsZ[1] = 10; // 10
ElmsZ[2] = 10; // 10
ElmsZ[3] = 5;  // 5

// Creation of the domain
// Points
For i In {1:4}
    Point(4*(i-1)+1) = { X[1], Y[i], Z[1], lc  }; 
    Point(4*(i-1)+2) = { X[2], Y[i], Z[1], lc  };
    Point(4*(i-1)+3) = { X[3], Y[i], Z[1], lc  };
    Point(4*(i-1)+4) = { X[4], Y[i], Z[1], lc  };
EndFor

// Lines
// horizontal lines
For i In {1:4}
    Line(3*(i-1)+1) = {4*(i-1)+1, 4*(i-1)+2};
    Line(3*(i-1)+2) = {4*(i-1)+2, 4*(i-1)+3};
    Line(3*(i-1)+3) = {4*(i-1)+3, 4*(i-1)+4};
EndFor

// vertical lines
For i In {1:4}
    Line(3*(i-1)+13) = {i,i+4};
    Line(3*(i-1)+14) = {i+4, i+8};
    Line(3*(i-1)+15) = {i+8, i+12};
EndFor

// Delimitation of surfaces
Curve Loop(1) = {1,16,-4,-13};
Curve Loop(2) = {2,19,-5,-16};
Curve Loop(3) = {3,22,-6,-19};

Curve Loop(4) = {4,17,-7,-14};
Curve Loop(5) = {5,20,-8,-17};
Curve Loop(6) = {6,23,-9,-20};

Curve Loop(7) = {7,18,-10,-15};
Curve Loop(8) = {8,21,-11,-18};
Curve Loop(9) = {9,24,-12,-21};

// Mesh subdivision
// othogonal to X axis
Transfinite Line{1,4,7,10} = ElmsX[1]+1 ;
Transfinite Line{2,5,8,11} = ElmsX[2]+1 ;
Transfinite Line{3,6,9,12} = ElmsX[3]+1 ;
// othogonal to Y axis
Transfinite Line{13,16,19,22} = ElmsY[1]+1 ;
Transfinite Line{14,17,20,23} = ElmsY[2]+1 ;
Transfinite Line{15,18,21,24} = ElmsY[3]+1 ;

// Creating Surfaces
For s In {1:9}
    Plane Surface(s) = {s};
    Transfinite Surface{s};
    Recombine Surface{s};
EndFor

// Extruding surfaces
// soil below the column (layer 1-2-3)
Extrude{0,0,Z(1)} {Surface{1:9}; Layers{ElmsZ[1]}; Recombine;}
Extrude{0,0,Z(1)} {Surface{222,200,178,156,134,112,90,68,46}; Layers{ElmsZ[1]}; Recombine;}
Extrude{0,0,Z(1)} {Surface{244,266,288,310,332,354,376,398,420}; Layers{ElmsZ[1]}; Recombine;}

// soil around column + column layer 4
Extrude{0,0,Z(2)} {Surface{618,596,574,552,530,508,486,464,442}; Layers{ElmsZ[2]}; Recombine;}

// soil around column + column layer 5
Extrude{0,0,Z(2)} {Surface{640,662,684,706,728,750,772,794,816}; Layers{ElmsZ[2]}; Recombine;}

// soil around column + column layer 6
Extrude{0,0,Z(2)} {Surface{838,860,882,904,926,948,970,992,1014}; Layers{ElmsZ[2]}; Recombine;}

// column above the ground
Extrude{0,0,Z(3)} {Surface{1124}; Layers{ElmsZ[3]}; Recombine;}

// Physical volumes
// physical volume below the column
Physical Volume("Soil_layer_1", 1) = {1:9};
Physical Volume("Soil_layer_2", 2) = {10:18};
Physical Volume("Soil_layer_3", 3) = {19:27};

Physical Volume("Soil_layer_4", 4) = {28:31,33:36};
Physical Volume("Soil_layer_5", 5) = {37:40,42:45};
Physical Volume("Soil_layer_6", 6) = {46:49,51:54};

// physical volume of the column
Physical Volume("building",7) = {32,41,50,55};

// Physical surfaces
// X0
Physical Surface("X0",1) = {45,111,177,419,353,287,617,551,485,639,705,771,837,903,969,1035,1101,1167};
// Xm
Physical Surface("Xm",2) = {81,147,213,235,301,367,433,499,565,741,807,675,1005,939,873,1071,1137,1203};
// Y0
Physical Surface("Y0",3) = {33,55,77,407,385,363,605,583,561,627,649,671,825,847,869,1023,1045,1067};
// Ym
Physical Surface("Ym",4) = {173,195,217,283,261,239,481,459,437,767,789,811,965,987,1009,1163,1185,1207};
// Z0
Physical Surface("Z0",5) = {1,2,3,4,5,6,7,8,9};
// Zm
Physical Surface("Zm",6) = {1036,1058,1080,1102,1146,1168,1190,1212};

// X0b
Physical Surface("X0b",7) = {697,895,1093,1233};
// Xmb
Physical Surface("Xmb",8) = {719,917,1115,1225};
// Y0b
Physical Surface("Y0b",9) = {657,855,1053,1221};
// Ymb
Physical Surface("Ymb",10) = {723,921,1119,1229};
// Z0b
Physical Surface("Z0b",11) = {530};
// Zmb
Physical Surface("Zmb",12) = {1234};

Mesh 3 ;
Save "Pillar.msh" ;