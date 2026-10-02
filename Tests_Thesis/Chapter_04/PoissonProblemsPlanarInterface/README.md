# Poisson problems: planar interface

Thesis section: 4.4.2.

Studies Poisson convergence on [0,2] x [0,1] x [0,1], divided by a planar nonconforming interface. The manufactured solution u = cos(pi*y) cos(pi*z) (2*x - x^2 + sin(pi*x)) is used to measure broken L2 and H1 errors for linear and quadratic hexahedra and compare segment-based, element-based, and RBF quadrature.

Run `main` after `initGReS`. Case-specific inputs and helpers are in `Input/` and `Utils/`. Results are written to `Output/`, whose previous contents are replaced on each run.
