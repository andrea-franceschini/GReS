// Gmsh project modified to refine the first line of elements above the bottom boundary

Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

nTip = 6;
Xsize = 1;
Ysize = 1;
nX = 61;
nY = 31;

tipSize = nTip/(nX-1);

// Define points
Point(1) = {0, 1, 0};
Point(2) = {tipSize, 1, 0};
Point(3) = {Xsize-tipSize, 1, 0};
Point(4) = {Xsize, 1, 0};

Point(5) = {Xsize, 1+Ysize, 0};
Point(6) = {Xsize-tipSize, 1+Ysize, 0};
Point(7) = {tipSize, 1+Ysize, 0};
Point(8) = {0, 1+Ysize, 0};


// Define lines
Line(1) = {1,2}; 
Line(2) = {2,3}; 
Line(3) = {3,4};
Line(4) = {4,5}; 
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};
Line(9) = {2,7};
Line(10) = {3,6};

// Define curve loop and surface
Curve Loop(1) = {1,9,7,8};
Curve Loop(2) = {2,10,6,-9};
Curve Loop(3) = {3,4,5,-10};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// Apply transfinite meshing
Transfinite Line{1} = nTip+1;
Transfinite Line{2} = nX;
Transfinite Line{3} = nTip+1;
Transfinite Line{4} = nY;
Transfinite Line{5} = nTip+1;
Transfinite Line{6} = nX;
Transfinite Line{7} = nTip+1;
Transfinite Line{8} = nY;
Transfinite Line{9} = nY;
Transfinite Line{10} = nY;

//Transfinite Surface {1,2,3};
//Recombine Surface {1,2,3};

// Define physical groups
Physical Curve("Interface_top2bot", 1) = {1,2,3};
Physical Curve("Top_load", 2) = {5,6,7};
Physical Curve("Lateral_load", 3) = {8};
Physical Curve("Lateral_fixed", 4) = {4};
Physical Surface("Domain_1", 1) = {1,2,3};

// Mesh and save
Mesh 2;

