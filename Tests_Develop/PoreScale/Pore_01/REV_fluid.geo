Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
SetFactory("OpenCASCADE");

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:

// create box 
Box(1) = {0,0,0,1,1,1};
Sphere(2) = {0,0,0,0.7};
Sphere(3) = {1,1,1,0.5};
BooleanDifference{ Volume{1}; Delete;}{ Volume{2}; Delete;};
BooleanDifference{ Volume{1}; Delete;}{ Volume{3}; Delete;};

//Recombine Surface {:};
MeshSize {:} = 0.1;
Mesh 3;
Save "Fluid.msh";
