Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

nX = 4;
nY = 21; // DEVE ESSERE DISPARI PER AVERE IL NODE IN Y=5
nZ = 30;

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:

X = 2.5;
Y = 10;
Z = 15;
xmin = 2.5;

Point(1) = {xmin, 0, 0, lc};
Point(2) = {xmin+X, 0, 0, lc};
Point(3) = {xmin+X, Y, 0, lc};
Point(4) = {xmin, Y, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line{1} = nX;
Transfinite Line{2} = nY;
Transfinite Line{3} = nX;
Transfinite Line{4} = nY;

Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Transfinite Surface {1}; // structured grid
Recombine Surface {1}; // using hexahedra

// extruding mesh along the z direction
Extrude {0, 0, 15} { Surface{1}; Layers{nZ}; Recombine;}
Physical Volume("Right", 1) = {1};
Physical Surface("Right_axial",2) = {17};
Physical Surface("Right_shear",3) = {26};
Physical Surface("Right_bottom",4) = {1};
Physical Surface("Right_interf",1) = {25};


Mesh 3;
Save "RightBlock_hexa.msh" ;
