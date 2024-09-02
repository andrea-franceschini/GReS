Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
SetFactory("OpenCASCADE");

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:


// create box 
//Box(1) = {0,0,-125,1000,1000,0};
Box(1) = {0,0,-200,1000,1000,-100};
//Physical Volume("Overburden", 1) = {1};
Physical Volume("Underburden",1)={1};
Physical Surface("Int_bottom",1) = {6};
Physical Surface("Bottom",2) = {5};

//Recombine Surface {:};
MeshSize {:} = 80;
Mesh 3;
Save "Underburden.msh";
