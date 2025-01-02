Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
SetFactory("OpenCASCADE");

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:


// create box 
Box(1) = {0,0,0,3,3,0.5};
Physical Volume("cube", 1) = {1};
Physical Surface("bottom",1) = {5};
Physical Surface("top",2) = {6};
//Physical Surface("lateral_normY",3) = {3,4};
Transfinite Surface {1,2,3,4,5,6};


//Recombine Surface {:};
MeshSize {:} = 0.2;
Mesh 3;
Save "Cube_hexa.msh";
