Mesh.Format = 1;
Mesh.MshFileVersion = 2.2;

// base square
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {0, 1, 0, 1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// 1 element in x and y
Transfinite Line{1,2,3,4} = 2;

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Surface{1}; // structured triangles

// extrusion in z with 10 layers
Extrude {0, 0, 10} {
  Surface{1};
  Layers{10};
}

// physical groups
Physical Volume("column", 1) = {1};
Physical Surface("bot",1) = {1};
Physical Surface("top",2) = {26};
Physical Surface("south",3) = {13};
Physical Surface("north",4) = {21};
Physical Surface("west",5) = {25};
Physical Surface("east",6) = {17};

Mesh 3;
Save "Column_tetra.msh";