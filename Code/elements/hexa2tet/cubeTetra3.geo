Mesh.Format = 1;
Mesh.MshFileVersion = 2.2; 

h = 20;

Point(1) = {0,0,0,h};
Point(2) = {0.5,0,0,h};
Point(3) = {0.8,0.7,0,h};
Point(4) = {0.3,0.4,0,h};

Line(1) = {1, 2};
Line(2) = {2,3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};
Transfinite Surface {1};
//Recombine Surface {1};

Extrude {0.5,0.3,0.6} {Surface{1}; Layers{1};}

Physical Volume("Domain_Left", 1) = {1};
Physical Surface("load",1) = {26};


Mesh 3;
Save "Tetra_3.msh";

