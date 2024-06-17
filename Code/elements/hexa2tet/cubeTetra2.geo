Mesh.Format = 1;
Mesh.MshFileVersion = 2.2; 

s = 10;
h = s+1;
Point(1) = {0,0,0,h};
Point(2) = {s,0,0,h};
Point(3) = {s,1,0,h};
Point(4) = {0,1,0,h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};
Transfinite Surface {1};
//Recombine Surface {1};

Extrude {0,0,1} {Surface{1}; Layers{1};}

Physical Volume("Domain_Left", 1) = {1};
Physical Surface("load",1) = {26};


Mesh 3;
Save "Tetra_5.msh";
