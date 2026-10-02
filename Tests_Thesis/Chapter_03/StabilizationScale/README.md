# Stabilized piecewise constant multiplier

Thesis section: 3.2.3.2.

Solves a manufactured Poisson problem on two stacked subdomains of the unit cube, with u = sin(pi*x) sin(pi*y) (1/2 + z). A sweep of the P0 traction-jump stabilization scale examines solvability, multiplier accuracy, and the excessive smoothing caused by large stabilization parameters.

Run `main` after `initGReS`. Case-specific inputs and helpers are in `Input/` and `Utils/`. Results are written to `Output/`, whose previous contents are replaced on each run.
