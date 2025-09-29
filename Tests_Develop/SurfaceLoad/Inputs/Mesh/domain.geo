Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:
lc = 1e-2;

// Models dimension
rad_i = 40.;
rad_e = 250.;
lay(1) = 15;
lay(2) = 15;
lay(3) = 10;
lay(4) = 10;

// Mesh partition
elm_prop_i =  5;
elm_prop_e =  3*elm_prop_i;
elm_prop_z = 5;

// Model creation
// Points
frac = Sqrt(2)/2;
Point(1) = {         0,          0, 0, lc};
Point(2) = {     rad_i,          0, 0, lc};
Point(3) = {         0,      rad_i, 0, lc};
Point(4) = {frac*rad_i, frac*rad_i, 0, lc};

Point(5) = {     rad_e,          0, 0, lc};
Point(6) = {         0,      rad_e, 0, lc};
Point(7) = {frac*rad_e, frac*rad_e, 0, lc};

// Lines
Line(1) = {1, 2};
Circle(2) = {4,1,2};
Circle(3) = {4,1,3};
Line(4) = {1, 3};
Transfinite Curve {1,4} = Round(rad_i/elm_prop_i)+1 Using Progression 1;
Transfinite Curve {2,3} = Round(((Pi*rad_i)/4)/elm_prop_i)+1 Using Progression 1;

Line(5) = {2, 5};
Circle(6) = {7,1,5};
Circle(7) = {7,1,6};
Line(8) = {3,6};
Transfinite Curve {5,8} = Round((rad_e-rad_i)/elm_prop_e)+1 Using Progression 1.1;
Transfinite Curve {6,7} = Round(((Pi*(rad_e-rad_i))/4)/elm_prop_e)+1 Using Progression 1.1;

Curve Loop(1) = {1,-2,3,-4};
Curve Loop(2) = {5,-6,7,-8,-3,2};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Recombine Surface {1,2}; // using hexahedra

// extruding mesh along the z direction
Extrude {0, 0, lay(1)} { Surface{1,2}; Layers{Round(lay(1)/elm_prop_z)+1}; Recombine;}
Extrude {0, 0, lay(2)} { Surface{30,62}; Layers{Round(lay(2)/elm_prop_z)+1}; Recombine;}
Extrude {0, 0, lay(3)} { Surface{84,116}; Layers{Round(lay(3)/elm_prop_z)+1}; Recombine;}
Extrude {0, 0, lay(4)} { Surface{138,170}; Layers{Round(lay(4)/elm_prop_z)+1}; Recombine;}

Physical Volume("layer1", 1) = {1,2};
Physical Volume("layer2", 2) = {3,4};
Physical Volume("layer3", 3) = {5,6};
Physical Volume("layer4", 4) = {7,8};
Physical Surface("topI",5) = {192};
Physical Surface("topE",6) = {224};
Physical Surface("bot",7) = {1,2};
Physical Surface("latX",8) = {17,71,125,179,41,95,149,203};
Physical Surface("latY",9) = {29,83,137,191,53,107,161,215};
Physical Surface("latC",10) = {45,99,153,207,49,103,157,211};

Mesh 3;
Save "QrCircle.msh";