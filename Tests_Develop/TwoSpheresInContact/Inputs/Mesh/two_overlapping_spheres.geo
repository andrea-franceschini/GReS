SetFactory("OpenCASCADE");

// -------------------------
R     = 1.0;
lc    = 0.3;
alpha = 0.20;        // overlap percentage of diameter

D  = 2*R;
h  = alpha * D;      // desired overlap thickness

// distance between centers = 2R - h
d = 2*R - h;

// Place spheres so they overlap
Sphere(1) = {-d/2, 0, 0, R};   // Left sphere
Sphere(2) = { d/2, 0, 0, R};   // Right sphere

left_box_x0 = -d/2 + R*(1 - alpha);
right_box_x0 = d/2 - R*(1 - alpha);

// Box(x) = {x0, y0, z0, dx, dy, dz}
Box(3) = {left_box_x0, -2*R, -2*R, 2*R, 4*R, 4*R};  // Box to subtract from Left sphere
Box(4) = {right_box_x0, 2*R, 2*R, -2*R, -4*R, -4*R}; // Box to subtract from right sphere

// Subtract box from right sphere
left_squished[] = BooleanDifference{ Volume{1}; Delete; }{ Volume{3}; Delete; }; // Here, Delete; ensures that the original volumes are deleted
right_squished[] = BooleanDifference{ Volume{2}; Delete; }{ Volume{4}; Delete; };

// -------------------------
// Make interface conforming
// -------------------------
volumes[] = BooleanFragments{ Volume{1, left_squished[0]}; Delete; }{ Volume{2, right_squished[]}; Delete; };

// -------------------------
Physical Volume("LeftSphere")  = {volumes[0]};
Physical Volume("RightSphere") = {volumes[1]};

// -------------------------
left_surfaces[]  = Boundary{ Volume{volIDs[0]}; };
right_surfaces[] = Boundary{ Volume{volIDs[1]}; };

Physical Surface("LeftSphereSurface")  = {left_surfaces[]};
Physical Surface("RightSphereSurface") = {right_surfaces[]};

// -------------------------
Mesh.CharacteristicLengthMax = lc;
Mesh 3;
Save "two_overlapping_spheres.vtk";