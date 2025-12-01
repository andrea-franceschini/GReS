clear
clc
testFiles = {fullfile('Mesh','testMesh.m');...
             fullfile('Simparam','testSimparam.m');...
             fullfile('Materials','testMaterials.m')
             fullfile('BoundaryConditions','testBoundaries.m')
             };

results = runtests(testFiles);

if any([results.Failed])
  error("Some test not passed");
else
  disp("All test run successfully")
end