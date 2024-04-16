// Gmsh project created on Fri Mar 15 13:58:35 2024

Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

lc = 0.2;
// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:

Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};

// Circle center
Point(3) = {0.5, 0, -0.4, lc};
Circle(1) = {1, 3, 2};

Extrude {0, 1, 0} {Line{1}; Layers{5}; Recombine;}

//Plane Surface(1) = {1};
//Transfinite Surface {1};

Physical Surface("Domain_test",1)={1};

Mesh 2;
Save "Mesh_coarse_HEXA.msh";

