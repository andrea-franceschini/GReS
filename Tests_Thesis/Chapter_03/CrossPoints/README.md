# Cross points

Thesis section: 3.3.3.

Partitions the unit cube into eight blocks with alternating 3-by-3-by-3 and 4-by-4-by-4 meshes. The exact field u = z gives unit flux on horizontal interfaces and zero flux elsewhere. The test compares nodal and P0 multiplier constraints at cross-points and cross-lines, including configurations with seven or twelve active interfaces.

Run `main` after `initGReS`. Case-specific inputs and helpers are in `Input/` and `Utils/`. Results are written to `Output/`, whose previous contents are replaced on each run.
