# Test case 1: Non uniform sliding

Thesis section: 5.5.4.1.

Validates embedded-fracture mechanics using a prism with an inclined crack that terminates inside the domain. Imposed vertical displacement generates a nonuniform sliding distribution. Displacement contours and fracture response provide a qualitative comparison with the benchmark of Borja and subsequent EFEM studies.

Run `main` after `initGReS`. Case-specific inputs and helpers are in `Input/` and `Utils/`. Results are written to `Output/`, whose previous contents are replaced on each run.
