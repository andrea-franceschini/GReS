Mesh.Format = 1;
Mesh.MshFileVersion = 2.2; 

<<<<<<< HEAD
h = 0.05;
=======
h = 0.25;
>>>>>>> origin/feature/cross_grid

Point(1) = {0,0,0,h};
Point(2) = {1,0,0,h};
Point(3) = {1,1,0,h};
Point(4) = {0,1,0,h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};
Physical Surface("domain",1) = {1};
<<<<<<< HEAD

Mesh 2;
Save "mesh2.msh";
=======
Recombine Surface{1};

Mesh 2;
Save "mesh2_quad.msh";
>>>>>>> origin/feature/cross_grid

