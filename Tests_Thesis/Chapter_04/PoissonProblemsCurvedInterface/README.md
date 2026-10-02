# Poisson problems: curved interface

Thesis section: 4.4.2.

Repeats the manufactured three-dimensional Poisson convergence study with a curved interface between nonconforming hexahedral meshes. Broken L2 and H1 errors assess the accuracy of segment-based, element-based, and RBF interface integration, including the effect of quadrature order on quadratic elements.

Run `main` after `initGReS`. Case-specific inputs and helpers are in `Input/` and `Utils/`. Results are written to `Output/`, whose previous contents are replaced on each run.
