// Gmsh project modified to refine the first line of elements above the bottom boundary

Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

Xsize = 1;
Ysize = 2;

nX = 200; // number of nodes along X direction
nY = 50; // number of nodes long Y direction
thin = 0.075; // size of first line of element above bottom boundary
dY = thin*Ysize/(nY-1);
// Define points
Point(1) = {0, 0, 0};
Point(2) = {Xsize, 0, 0};
Point(3) = {Xsize, Ysize/2, 0};
Point(4) = {Xsize, Ysize/2+dY, 0};
Point(5) = {Xsize, Ysize, 0};
Point(6) = {0,Ysize,0};
Point(7) = {0,Ysize/2+dY,0};
Point(8) = {0,Ysize/2,0};



// Define lines
Line(1) = {1, 2}; // Bottom boundary
Line(2) = {2, 3}; // Right thin line
Line(3) = {3, 8}; // i
Line(4) = {8, 1}; // Left boundary
Line(5) = {3,4};
Line(6) = {4,7};
Line(7) = {7,8};
Line(8) = {4,5};
Line(9) = {5,6};
Line(10) = {6,7};




// Define curve loop and surface
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {-3,5,6,7};
Curve Loop(3) = {-6,8,9,10};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// Apply transfinite meshing
Transfinite Line{1,3,6,9} = nX;
Transfinite Line{2,4,8,10} = nY;
Transfinite Line{5,7} = 2;
Transfinite Surface {1,2,3};

Recombine Surface {1,2,3};

// Define physical groups
Physical Curve("Bottom_fixed", 1) = {1};
Physical Curve("Lateral_fixed", 2) = {2,5,8};
Physical Surface("Domain_1", 1) = {1,3};
Physical Surface("Thin_layer",2) = {2};
// Mesh and save
Mesh 2;
Save "Block_hexa.msh";

