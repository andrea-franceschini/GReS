// Gmsh project created on Fri Mar 15 13:58:35 2024

Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

NX = 32;
NY = 32;
xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;

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

Transfinite Line{1} = NX+1;
Transfinite Line{2} = NY+1;
Transfinite Line{3} = NX+1;
Transfinite Line{4} = NY+1;
Transfinite Surface {1};
Recombine Surface {1};

Physical Curve("Interface",1) = {3};
Physical Curve("Fix_x",2) = {2,4};
Physical Curve("Fix_y",3) = {1};
Physical Surface("slave_bot",1)={1};

//Mesh 2;
//Save "slave.msh";

