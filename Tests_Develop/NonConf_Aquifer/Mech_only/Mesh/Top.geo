Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
SetFactory("OpenCASCADE");

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:


// create box 
Box(1) = {0,0,0,500,500,-200};
//Box(2) = {0,0,-300,1000,1000,-175};
Physical Volume("Top", 1) = {1};
//Physical Volume("Underburden",2)={2};
Physical Surface("Int_top",1) = {5};
Physical Surface("Top_load",2) = {6};


//Recombine Surface {:};
MeshSize {:} = 400;
Mesh 3;
Save "Top.msh";
