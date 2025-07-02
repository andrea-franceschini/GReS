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
Sphere(2) = {0.1,0.1,0,0.5};
Sphere(3) = {1,1,1,0.3};
Sphere(4) = {0.85,0.4,0.5,0.3};
Sphere(5) = {0,1,0.9,0.3};
BooleanIntersection(6) = { Volume{1};}{ Volume{2}; Delete;};
BooleanIntersection(7) = { Volume{1};}{ Volume{3}; Delete;};
BooleanIntersection(8) = { Volume{1};}{ Volume{4}; Delete;};
BooleanIntersection(9) = { Volume{1}; Delete;}{ Volume{5}; Delete;};

//Recombine Surface {:};
MeshSize {:} = 0.3;
Mesh 3;
Save "Grain_3.msh";
