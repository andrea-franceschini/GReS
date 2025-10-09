// 3d_sphere_fast.geo — quick sphere mesh with modest center refinement
Mesh.Format = 1;
Mesh.MshFileVersion = 2.2;

SetFactory("OpenCASCADE");

// basic sizing (coarser for speed)
R = 1.0;
s_min = 0.05;   // smallest element size near center (coarser -> faster)
s_max = 0.01;    // largest element size near boundary (coarser -> faster)
alpha = 2.0;     // controls concentration near center

// fallback global sizes
Mesh.CharacteristicLengthMax = 0.1;
Mesh.CharacteristicLengthMin = 0.03;

// geometry
Sphere(1) = {0, 0, 0, R};
Physical Volume("SphereVolume") = {1};
Physical Surface("SphereSurface") = {1};

// simple radial MathEval field: size = s_min + (s_max - s_min)*(r/R)^alpha
Field[1] = MathEval;
Field[1].F = sprintf("%g + (%g - %g)*pow(sqrt(x*x+y*y+z*z)/%g, %g)", s_min, s_max, s_min, R, alpha);
Background Field = 1;

// a recommended (sometimes faster) 3D algorithm: uncomment if you want to try it
// Mesh.Algorithm3D = 1; // 1 = Delaunay, 4 = Frontal (default on many installs)

// mesh and save
Mesh 3;
Save "3d_sphere.msh";



