# Test case 2: Single crack under compression

Thesis section: 5.4.3.2.

A finite crack inclined by 20 degrees is subjected to uniform uniaxial compression under plane-strain conditions. With fracture length 2, friction angle 30 degrees, and zero cohesion, numerical normal traction and tangential gap profiles are compared with analytical expressions. Mortar tying outside the crack represents the surrounding intact material.

Run `main` after `initGReS`. Case-specific inputs and helpers are in `Input/` and `Utils/`. Results are written to `Output/`, whose previous contents are replaced on each run. This case uses its full configuration in normal runs and is skipped entirely by the suite when `Smoke=true`.
