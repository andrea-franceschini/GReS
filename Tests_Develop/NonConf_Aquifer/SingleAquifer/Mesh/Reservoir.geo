Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
//SetFactory("OpenCASCADE");

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:


// create box
// Different boxes are not glued toghether! 
lc = 100;
Point(1) = {0,0,-200,lc};
Point(2) = {1000,0,-200,lc};
Point(3) = {1000,1000,-200,lc};
Point(4) = {0,1000,-200,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Extrude{0,0,40} {Surface{1};Layers{2};}
Extrude{0,0,20} {Surface{26};Layers{3};}
Extrude{0,0,40} {Surface{48};Layers{2};}


Physical Volume("Bot", 3) = {1};
Physical Volume("Mid", 2) = {2};
Physical Volume("Top", 1) = {3};
Physical Surface("Res_top",1) = {70};
Physical Surface("Res_bot",2) = {1};


//Recombine Surface {:};
//MeshSize {:} = 50;
Mesh 3;
Save "Reservoir.msh";
