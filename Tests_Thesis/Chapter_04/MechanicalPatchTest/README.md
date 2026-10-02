# Mechanical patch test (supplementary)

Thesis section: 4.4.3 (supplementary mechanical variant).

Applies uniform compression to two vertically stacked elastic blocks with nonconforming meshes. Zero Poisson ratio gives an exact uniaxial constant-strain solution. Displacement errors compare segment-based, element-based, and RBF mortar integration; this supplements the scalar Poisson patch test described in the thesis.

Run `main` after `initGReS`. Case-specific inputs and helpers are in `Input/` and `Utils/`. Results are written to `Output/`, whose previous contents are replaced on each run.
