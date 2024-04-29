Mesh.Format = 1;
Mesh.MshFileVersion = 2.2; 

h = 0.5;

Point(1) = {1,0,0,h};
Point(2) = {2,0,0,h};
Point(3) = {2,1,0,h};
Point(4) = {1,1,0,h};
Point(5) = {0.75,0.5,0,h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Circle(4) = {4,5,1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};
//Transfinite Surface {1};
Recombine Surface {1};

Extrude {0,0,1} {Surface{1}; Layers{12}; Recombine;}

Physical Volume("Domain_Right", 1) = {1};
Physical Surface("Interface_Right",1) = {25};
Physical Surface("Dirichlet",2) = {1,13,17,21,26};


Mesh 3;
Save "RightBlock_curve.msh";

