clear
clc

v = gresLog().getVerbosity;
gresLog().setVerbosity(-2);
testFiles = {fullfile('Mesh','testMesh.m');...
             fullfile('Simparam','testSimparam.m');...
             fullfile('Materials','testMaterials.m')
             fullfile('BoundaryConditions','testBoundaries.m')
             fullfile('PatchTestMechanics','testShearPatch.m')
             };

results = runtests(testFiles);
gresLog().setVerbosity(v);
if any([results.Failed])
  error("Some test not passed");
else
  disp("All test passed")
end

