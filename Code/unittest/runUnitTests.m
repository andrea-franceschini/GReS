clear
clc


testFiles = {fullfile('Mesh','testMesh.m');...
             fullfile('DoFManager','testDoFManager.m');...
             fullfile('Simparam','testSimparam.m')};

results = runtests(testFiles);

if any([results.Failed])
  error("Some unit tests not passed");
else
  disp("All test run successfully")
end
