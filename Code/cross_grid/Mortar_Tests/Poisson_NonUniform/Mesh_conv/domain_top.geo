Mesh.Format = 1; // msh output format
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use

// Mesh size parameters
RAT = 3;
// Number of divisions along x and y
NX = 40;
NY = Floor(NX/2);

// Compute progression factor to smoothly transition element sizes along x
ProgressionX = (RAT)^(1.0 / (NX - 1)); // Geometric progression from NL to NR


// Geometry definition
Point(1) = {0, 0.5, 0};
Point(2) = {1, 0.5, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {4, 3};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, -3, 4};
Plane Surface(1) = {1};

// Apply transfinite lines with non-uniform distribution
Transfinite Line{1,3} = NX+1 Using Progression ProgressionX;  // Non-uniform x-mesh
Transfinite Line{2,4} = NY+1;  // Uniform mesh along y

// Keep the structured mesh
Transfinite Surface{1};


Physical Curve("Interface_top",1) = {1};
Physical Curve("External_boundary",2) = {2,3,4};
Physical Surface("Domain_2",1)={1};


