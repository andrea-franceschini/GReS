% % % % clear
% % % % clc
% % % % 
% % % % v = gresLog().getVerbosity;
% % % % gresLog().setVerbosity(-2);
% % % % testFiles = {fullfile('Mesh','testMesh.m');...
% % % %              fullfile('Simparam','testSimparam.m');...
% % % %              fullfile('Materials','testMaterials.m')
% % % %              fullfile('BoundaryConditions','testBoundaries.m')
% % % %              fullfile('PatchTestMechanics','testShearPatch.m')
% % % %              };
% % % % 
% % % % results = runtests(testFiles);
% % % % gresLog().setVerbosity(v);
% % % % if any([results.Failed])
% % % %   error("Some test not passed");
% % % % else
% % % %   disp("All test passed")
% % % % end
% % % % 


% runUnitTests.m
%
% Run GReS unit tests in parallel.

testPath = fullfile(gres_root(),'Code','unittest');

pool = gcp('nocreate');
if isempty(pool)
    parpool('Processes', 2);
end

results = runtests(fullfile(testPath, 'UnitTests.m'),'UseParallel', true);
disp(results);
if any([results.Failed])
  error('GReS:UnitTestsFailed','One or more unit tests failed.');
else
  disp('All unit tests passed.');
end