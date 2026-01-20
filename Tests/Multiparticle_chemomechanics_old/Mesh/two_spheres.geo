// ----------------------------------------------------------------------------
// 3D mesh of two tangent spheres (each tagged separately)
// ----------------------------------------------------------------------------
SetFactory("OpenCASCADE");

// mesh sizing
R  = 1.0;
lc = 0.2;   // increase for coarser, decrease for finer mesh

// define the two spheres
Sphere(1) = {-1, 0, 0, R};    // left sphere centered at (-1, 0, 0)
Sphere(2) = { 1.000001, 0, 0, R};    // right sphere centered at (1, 0, 0)

// make sure geometry is explicit (keeps separate volumes even if tangent)
BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// ---------------------------------------------------------------------------
// Tagging volumes and surfaces individually
// ---------------------------------------------------------------------------

// list of all volumes and surfaces to verify
// (you can run: gmsh -0 two_spheres_boolean.geo to see their IDs)
volIDs[] = Volume{:};
surfIDs[] = Surface{:};

// For safety, we’ll assume the first two volumes correspond to the two spheres.
// If BooleanFragments creates more than two volumes (unlikely here), you can check
// with gmsh -0 and adjust the indices below accordingly.
Physical Volume("SphereLeft")  = {volIDs[0]};
Physical Volume("SphereRight") = {volIDs[1]};

// The outer surfaces are usually 1 and 2 (inner fragments possible if Boolean adds interface).
// We’ll conservatively assign surfaces belonging to each sphere.
left_surfaces[]  = Boundary{ Volume{volIDs[0]}; };
right_surfaces[] = Boundary{ Volume{volIDs[1]}; };

Physical Surface("LeftSphereSurface")  = {left_surfaces[]};
Physical Surface("RightSphereSurface") = {right_surfaces[]};

// optional: if tangency creates a shared vertex
// you can find its id with gmsh -0 and then add:
// Physical Point("ContactPoint") = {vertex_id};

// ---------------------------------------------------------------------------
// mesh control
// ---------------------------------------------------------------------------
Mesh.CharacteristicLengthMax = lc;
Mesh.CharacteristicLengthMin = lc/3;

// generate and save directly to VTK
Mesh 3;
Save "two_spheres.vtk";
