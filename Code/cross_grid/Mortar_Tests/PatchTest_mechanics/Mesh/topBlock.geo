// Gmsh project modified to refine the first line of elements above the bottom boundary

Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

Xsize = 1;
Ysize = 0.5;
nX = 16;
nY = 16;
//thin = 1; // size of first line of element above bottom boundary
//dY = thin*Ysize/(nY-1);
// Define points
Point(1) = {0, 0.5, 0};
Point(2) = {Xsize, 0.5, 0};
Point(3) = {Xsize, 0.5+Ysize, 0};
Point(4) = {0, 0.5+Ysize, 0};

// Define lines
Line(1) = {1, 2}; 
Line(2) = {2, 3}; 
Line(3) = {3, 4}; 
Line(4) = {4, 1}; 


// Define curve loop and surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Apply transfinite meshing
Transfinite Line{1,3} = nX;
Transfinite Line{2,4} = nY;
Transfinite Surface {1};
Recombine Surface {1};

// Define physical groups
Physical Curve("Interface_top2bot", 1) = {1};
Physical Curve("Top_load", 2) = {3};
Physical Curve("Lateral_load", 3) = {4};
Physical Curve("Lateral_fixed", 4) = {2};
Physical Surface("Domain_1", 1) = {1};

// Mesh and save
//Mesh 2;
//Save "top.msh";
