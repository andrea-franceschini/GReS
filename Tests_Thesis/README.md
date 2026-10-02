# Thesis numerical experiments

Standalone GReS version to reproduce the test case of Daniele Moretto's PhD thesis.

## Run

Initialize GReS from the root folder:

```matlab
initGReS
compileAll
```

Then move to the folder containing the test cases
```matlab
cd('currentFolder/Thesis_tests')
report = runAllTests('Smoke',true);
report = runAllTests();         % full suite, serial
```
The tests use the GReS checkout selected by `initGReS`. GReS MEX functions
must already be compiled for the current platform. Smoke mode reduces selected
mesh refinements in supported cases. With `Smoke=true`, the contact comparison,
single crack under compression, and graben-horst cases are skipped entirely:
the global runner selects 17 cases instead of 20. With `Smoke=false` (the
default), all 20 cases are selected. Full Q2 studies and the Sneddon 405-by-405 mesh can be costly.
The runnable cases generate their meshes in GReS without Python, Gmsh, or MRST.

## Organization

Every case has `main.m`, `Utils/`, `Input/`, and `Output/`. The experiment is in
`main.m`; `Utils` contains only additional functions specific to that case.
Mesh and boundary helpers used by a test case live in the `Utils` folder of that case; the suite-level `Utils` holds only the runner and the optional validation mesh helper.

```matlab
runAllTests('ListOnly',true)                  % inventory
runAllTests('Chapter',4,'Smoke',true)        % one chapter
runAllTests('Cases',"PatchTest")             % selected case
```
Each chapter also has a `runAllTests.m`. Run it from that chapter folder. Cases
run in sequence; one failure is recorded and the runner continues to the next.

Each case writes to its own `Output/` folder. Starting that case again deletes the
previous contents of that folder, then saves `run.log`, `status.mat`, and (on
success) `result.mat` alongside the files written by `OutState`. The suite
summary is overwritten in the top-level `Output/report.csv` and `report.mat`.
Input files stay in the case's `Input/` folder. 

Part of the reorganization of the folder has been carried out with the aid of AI.

