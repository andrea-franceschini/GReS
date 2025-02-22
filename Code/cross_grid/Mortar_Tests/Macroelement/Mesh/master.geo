// Gmsh project created on Fri Mar 15 13:58:35 2024

Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

N = 2;
xMin = 0;
xMax = 2;
yMin = 1;
yMax = 2;

Point(1) = {xMin, yMin, 0};
Point(2) = {xMax, yMin, 0};
Point(3) = {xMax, yMax, 0};
Point(4) = {xMin, yMax, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line{1} = N+1;
Transfinite Line{2} = 2;
Transfinite Line{3} = N+1;
Transfinite Line{4} = 2;
Transfinite Surface {1};
Recombine Surface {1};

Physical Curve("Interface",1) = {1};
Physical Curve("Fix",2) = {2,3,4};
Physical Surface("master_top",1)={1};

//Mesh 2;
//Save "slave.msh";


