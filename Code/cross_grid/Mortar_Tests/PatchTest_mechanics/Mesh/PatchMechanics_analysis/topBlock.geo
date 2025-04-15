// Gmsh project modified to refine the first line of elements above the bottom boundary

Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

Xsize = 1;
Ysize = 0.5;
nX = 80;
nY = 80;

// Define points
Point(1) = {0, Ysize, 0};
Point(2) = {0.5*Xsize, Ysize, 0};
Point(3) = {Xsize, Ysize, 0};
Point(4) = {Xsize, 2*Ysize, 0};
Point(5) = {0.5*Xsize, 2*Ysize, 0};
Point(6) = {0, 2*Ysize, 0};

// Define lines
Line(1) = {1, 2}; 
Line(2) = {2, 3}; 
Line(3) = {3, 4}; 
Line(4) = {4, 5}; 
Line(5) = {5, 6}; 
Line(6) = {6, 1}; 
Line(7) = {2, 5}; 

// Define curve loop and surface
Curve Loop(1) = {1, 7, 5, 6};
Plane Surface(1) = {1};

Curve Loop(2) = {2,3,4,-7};
Plane Surface(2) = {2};


// Apply transfinite meshing
Transfinite Line{1,2,4,5} = 0.5*nX+1;
Transfinite Line{3,6,7} = 0.5*nY+1;

Transfinite Surface {1};
Recombine Surface {1};
Transfinite Surface {2};
Recombine Surface {2};

// Define physical groups
Physical Curve("Interface_top2bot", 1) = {1,2};
Physical Curve("Top_load", 2) = {4,5};
Physical Curve("Lateral_load", 3) = {6};
Physical Curve("Lateral_fixed", 4) = {3};
Physical Surface("Domain_left", 1) = {1};
Physical Surface("Domain_right", 2) = {2};

// Mesh and save
//Mesh 2;
//Save "top.msh";
