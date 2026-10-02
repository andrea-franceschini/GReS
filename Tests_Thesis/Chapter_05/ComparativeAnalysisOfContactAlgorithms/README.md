# A comparative analysis of contact algorithms

Thesis section: 5.6; loading benchmark in 5.4.3.3 and 5.5.4.4.

Two elastic blocks separated by a frictional crack undergo compression and shear, including load reversal. Separate EFEM and mortar models compare contact algorithms under vertical and horizontal loading, before and after reversal. Fixed-parameter studies and augmentation sweeps report total nonlinear iterations and failed configurations.

Run `main` after `initGReS`. Both formulations run serially and write separate results under `Output/EFEM` and `Output/Mortar`. Normal runs use the full fixed-parameter study and augmentation sweep. The suite skips this case entirely when `Smoke=true`; the standalone entry has no smoke option.
