Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
// Mesh.Format = 16; // vtk output format

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:
lc = 1e-2;
filename = "dam.msh";
// filename = "dam.vtk";

// Geometry
//  XZ view
//                            top
//                      2-------------1
//                     /               \
//                    /      middle     \
//                   /       4---3       \
//                  /       /     \       \
//                 /       /       \       \
//     13---------5       7         8       6-----------14
//     |                 /           \                   |
//     |                /             \                  |
//     15-------9-----11---------------12------10-------16


// Models dimension
ang(1) = 30.;
ang(2) = 60.;

top(1) = 4; // 0.4 default base width
top(2) = 25;  // 1.0 default height
middle(1) = 2;  // 0.2 default top width
middle(2) = 23;  // 0.7 default middle height
lay(1) = 20;   // 0.3 default layer height
lay(2) = 10;   // 0.5 default layer thickness below base

angRad(1) = (90.-ang(1))*Pi/180.;
angRad(2) = (90.-ang(2))*Pi/180.;

sectLength = 7.5;  // 1.5 default section length in the y direction

// mesh dimensions
Elms[1] = 64;
Elms[2] = 32;

Elms[3] = 16;

Elms[4] = 32;
Elms[5] = 32;


Elms[6] = 16;


// Elms[1] = 32;
// Elms[2] = 16;

// Elms[3] = 8;

// Elms[4] = 16;
// Elms[5] = 16;


// Elms[6] = 8;

// Model creation
// Points
Point(1) = {    top(1)/2., 0.,    top(2), lc};
Point(2) = {   -top(1)/2., 0.,    top(2), lc};

Point(3) = { middle(1)/2., 0., middle(2), lc};
Point(4) = {-middle(1)/2., 0., middle(2), lc};

len = top(1)/2. + Tan(angRad(1))*(top(2)-lay(1));
Point(5) = { len, 0., lay(1), lc};
Point(6) = {-len, 0., lay(1), lc};

len = middle(1)/2. + Tan(angRad(2))*(middle(2)-lay(1));
Point(7) = { len, 0., lay(1), lc};
Point(8) = {-len, 0., lay(1), lc};

len = top(1)/2. + Tan(angRad(1))*top(2);
Point( 9) = { len, 0., 0., lc};
Point(10) = {-len, 0., 0., lc};

len = middle(1)/2. + Tan(angRad(2))*middle(2);
Point(11) = { len, 0., 0., lc};
Point(12) = {-len, 0., 0., lc};

len = Round(top(1)/2. + Tan(angRad(1))*top(2) + lay(2));
Point(13) = { len, 0., lay(1), lc};
Point(14) = {-len, 0., lay(1), lc};
Point(15) = { len, 0.,     0., lc};
Point(16) = {-len, 0.,     0., lc};

// Lines
Line(1) = { 15,  9};
Line(2) = {  9,  5};
Line(3) = {  5, 13};
Line(4) = { 13, 15};

Line( 5) = {  9, 11};
Line( 6) = { 11,  7};
Line( 7) = {  7,  5};

Line( 8) = { 11, 12};
Line( 9) = { 12,  8};
Line(10) = {  8,  7};

Line(11) = { 12, 10};
Line(12) = { 10,  6};
Line(13) = {  6,  8};

Line(14) = { 10, 16};
Line(15) = { 16, 14};
Line(16) = { 14,  6};

Line(17) = {  7,  3};
Line(18) = {  3,  1};
Line(19) = {  1,  5};

Line(20) = {  8,  4};
Line(21) = {  4,  3};

Line(22) = {  6,  2};
Line(23) = {  2,  4};

Line(24) = {  2,  1};

Transfinite Curve{2,4,6,9,12,15} = Elms[1]+1;
Transfinite Curve{17,19,20,22} = Elms[2]+1;

Transfinite Curve{8,10,21,24} = Elms[3]+1;

Transfinite Curve{5,7,11,13,18,23} = Elms[4]+1;
Transfinite Curve{-1,3,14,-16} = Elms[5]+1 Using Progression 1.1;

// Surfaces
Curve Loop(1) = {  8,  9, 10, -6};
Curve Loop(2) = {-10, 20, 21,-17};

Curve Loop(3) = {  1,  2,  3,  4};
Curve Loop(4) = { 14, 15, 16,-12};

Curve Loop(5) = {  5,  6,  7, -2};
Curve Loop(6) = { 11, 12, 13, -9};

Curve Loop(7) = { -7, 17, 18, 19};
Curve Loop(8) = {-13, 22, 23,-20};

Curve Loop(9) = {-21,-23, 24,-18};

For s In {1:9}
    Plane Surface(s) = {s};
    Transfinite Surface {s}; // structured grid
    Recombine Surface {s}; // using hexahedra
EndFor

// extruding mesh along the z direction
Extrude {0, sectLength, 0} { Surface{1:9}; Layers{Elms[6]}; Recombine;}

Physical Volume("mat1", 1) = {1,2};
Physical Volume("mat2", 2) = {3:9};

Physical Surface("latX0",1) = {103};
Physical Surface("latXM",2) = {89};

Physical Surface("latY0",3) = {1:9};
Physical Surface("latYM",4) = {46,68,90,112,134,156,178,200,222};

Physical Surface("top",5) = {85,107,177,191,217};
Physical Surface("bot",6) = {33,77,99,121,143};

Physical Surface("topA",7) = {85,177};
Physical Surface("topB",8) = {217};
Physical Surface("topC",9) = {107,191};

Mesh 3;
Save Sprintf(filename);