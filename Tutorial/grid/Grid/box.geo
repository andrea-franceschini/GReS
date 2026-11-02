Mesh.Format = 16; // vtk output format

Mesh.MshFileVersion = 2.2; // Required version for .msh file format

Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


Transfinite Line{1} = 2;
Transfinite Line{2} = 2;
Transfinite Line{3} = 2;
Transfinite Line{4} = 2;

Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Transfinite Surface {1}; // structured grid
Recombine Surface {1}; // using hexahedra

// extruding mesh along the z direction
Extrude {0, 0, 1} { Surface{1}; Layers{1}; Recombine;}

Physical Volume("box", 1) = {1};
Physical Surface("top",2) = {26};
Physical Surface("bottom",1) = {1};

Mesh 3;
Save "Box.vtk";
