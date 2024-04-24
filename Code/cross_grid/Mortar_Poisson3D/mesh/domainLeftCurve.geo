Mesh.Format = 1;
Mesh.MshFileVersion = 2.2; 

h = 0.1;

Point(1) = {0,0,0,h};
Point(2) = {1,0,0,h};
Point(3) = {1,0,1,h};
Point(4) = {0,0,1,h};
Point(5) = {1.3,0,0.5,h};
Line(1) = {1, 2};
Circle(2) = {2, 5, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};
//Transfinite Surface {1};
Recombine Surface {1};

Extrude {0,1,0} {Surface{1}; Layers{10}; Recombine;}

Physical Volume("Domain_Left", 1) = {1};
Physical Surface("Interface_Left",1) = {17};
Physical Surface("Dirichlet",2) = {1,13,25,21,26};


Mesh 3;
Save "domainLeftCurve_H1.msh";

