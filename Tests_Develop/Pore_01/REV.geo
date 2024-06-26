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
Sphere(2) = {0,0,0,0.2};
Sphere(3) = {1,1,1,0.4};
BooleanIntersection(4) = { Volume{1};}{ Volume{2}; Delete;};
BooleanIntersection(5) = { Volume{1}; Delete;}{ Volume{3}; Delete;};

//Mesh.RecombinationAlgorithm = 2;
//Transfinite Surface {:};
Recombine Surface {:};
//MeshSize {:} = 0.15;
//Mesh 3;
//RecombineMesh;
Mesh.SubdivisionAlgorithm = 3;
//RecombineMesh;
//Mesh.Algorithm = 8;
MeshSize {:} = 0.12;
Mesh 3;
