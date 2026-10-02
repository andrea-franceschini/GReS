# Computational cost

Thesis section: 4.3.5.

Compares analytical floating-point operation estimates for element-based Gauss-point projection and RBF mortar integration. The calculation varies the number of Gauss points and relative mesh refinement while accounting for projection iterations and RBF interpolation setup. The exported curves represent estimated cost, rather than measured solver timings.

Run `main` after `initGReS`. Case-specific inputs and helpers are in `Input/` and `Utils/`. Results are written to `Output/`, whose previous contents are replaced on each run.
