Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
SetFactory("OpenCASCADE");

// create box 
Box(1) = {0,0,0,1,1,5};
Physical Volume("cube") = {1};
Physical Surface("xmin") = {1};
Physical Surface("xmax") = {2};
Physical Surface("ymin") = {3};
Physical Surface("ymax") = {4};
Physical Surface("zmin") = {5};
Physical Surface("zmax") = {6};
Transfinite Surface {1,2,3,4,5,6};

MeshSize {:} = 0.2;
Mesh 3;
Save "column.msh";
