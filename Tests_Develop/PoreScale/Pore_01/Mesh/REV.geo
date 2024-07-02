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
BooleanIntersection(4) = { Volume{1};}{ Volume{2}; Delete;};
BooleanIntersection(5) = { Volume{1}; Delete;}{ Volume{3}; Delete;};
Physical Volume("g1", 1) = {4};
Physical Volume("g2",2)={5};
Physical Surface("fix",1) = {9,10,12,13,14,15,16,18};
Physical Surface("int1",2) = {11};
Physical Surface("int2",3) = {17};


//Recombine Surface {:};
MeshSize {:} = 40;
Mesh 3;
Save "Grain.msh";
