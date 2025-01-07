// Gmsh project modified to refine the first line of elements above the bottom boundary

Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

Xsize = 1;
Ysize = 1;
nX = 51;
nY = 26;
thin = 1; // size of first line of element above bottom boundary
dY = thin*Ysize/(nY-1);
// Define points
Point(1) = {0, 1, 0};
Point(2) = {Xsize, 1, 0};
Point(3) = {Xsize, 1+dY, 0};
Point(4) = {Xsize, 1+Ysize, 0};
Point(5) = {0, 1+Ysize, 0};
Point(6) = {0,1+dY,0};

// Define lines
Line(1) = {1, 2}; // Bottom boundary
Line(2) = {2, 3}; // Right thin line
Line(3) = {3, 4}; // i
Line(4) = {4, 5}; // Left boundary
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {3,6};

// Define curve loop and surface
Curve Loop(1) = {1, 2, 7, 6};
Curve Loop(2) = {-7,3,4,5};
Plane Surface(1) = {1};
Plane Surface(2) = {2};

// Apply transfinite meshing
Transfinite Line{1} = nX;
Transfinite Line{2} = 2;
Transfinite Line{3} = nY;
Transfinite Line{4} = nX;
Transfinite Line{5} = nY;
Transfinite Line{6} = 2;
Transfinite Line{7} = nX;
Transfinite Surface {1,2};
Recombine Surface {1,2};

// Define physical groups
Physical Curve("Interface_top2bot", 1) = {1};
Physical Curve("Top_load", 2) = {4};
Physical Curve("Lateral_load", 3) = {5};
Physical Curve("Lateral_fixed", 4) = {2,3};
Physical Surface("Domain_1", 1) = {1,2};

// Mesh and save
//Mesh 2;
//Save "top.msh";
