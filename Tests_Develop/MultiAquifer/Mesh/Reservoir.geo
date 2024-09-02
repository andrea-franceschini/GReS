Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
SetFactory("OpenCASCADE");

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:


// create box 
Box(1) = {0,0,-100,1000,1000,-25};
Box(2) = {0,0,-125,1000,1000,-50};
Box(3) = {0,0,-175,1000,1000,-25};
Physical Volume("L1", 1) = {1};
Physical Volume("L2", 2) = {2};
Physical Volume("L3", 3) = {3};
Physical Surface("Res_top",1) = {6};
Physical Surface("Res_bot",2) = {17};


//Recombine Surface {:};
MeshSize {:} = 20;
Mesh 3;
Save "Reservoir.msh";
