Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
SetFactory("OpenCASCADE");

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:

// create box 
Box(1) = {0,0,0,100,100,100};
Sphere(2) = {20,50,5,60};
Sphere(3) = {90,50,110,65};
BooleanDifference{ Volume{1}; Delete;}{ Volume{2}; Delete;};
BooleanDifference{ Volume{1}; Delete;}{ Volume{3}; Delete;};

Physical Volume("wat",1)={1};
Physical Surface("bottom",1) = {6};
Physical Surface("top",2) = {3};
Physical Surface("g1",3) = {5};
Physical Surface("g2",4) = {8};

//Recombine Surface {:};
MeshSize {:} = 8;
Mesh 3;
Save "Fluid.msh";
