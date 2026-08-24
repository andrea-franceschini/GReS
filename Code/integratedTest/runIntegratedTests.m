% % % % clear
% % % % clc
% % % % 
% % % % tWall = tic();      % start wall-clock timer
% % % % tCPU  = cputime;    % start CPU timer
% % % % 
% % % % v = gresLog().getVerbosity();
% % % % gresLog().setVerbosity(-2);
% % % % 
% % % % testFiles = {
% % % %     fullfile('Terzaghi','testTerzaghi.m')
% % % %     fullfile('SubDomains','testSubDomains.m')
% % % %     fullfile('Richards','testRichards.m')
% % % %     fullfile('CubeSedimentation','testSedimentation.m')
% % % %     fullfile('MortarConvergence','testMortarPoisson.m')
% % % %     fullfile('ConstantSliding','testConstantSliding.m')
% % % %     fullfile('SingleCrackCompressed','testSingleCrackCompressed.m')
% % % %     fullfile('ConstantSlidingEFEM','testConstantSlidingEFEM.m')
% % % % };
% % % % 
% % % % results = runtests(testFiles,'UseParallel', true);
% % % % 
% % % % elapsedWall = toc(tWall);
% % % % elapsedCPU  = cputime - tCPU;
% % % % 
% % % % fprintf("Elapsed wall-clock time: %1.2f s\n", elapsedWall);
% % % % fprintf("Elapsed CPU time:        %1.2f s\n", elapsedCPU);
% % % % 
% % % % gresLog().setVerbosity(v);
% % % % 
% % % % if any([results.Failed])
% % % %     error("Some tests did not pass");
% % % % else
% % % %     disp("All tests passed");
% % % % end


% runIntegratedTests.m
%
% Run GReS integrated tests in parallel.

testPath = fullfile( ...
    gres_root(), ...
    'Code', ...
    'integratedTest');

pool = gcp('nocreate');

if isempty(pool)
    parpool('Processes', 2);
end

results = runtests( ...
    fullfile(testPath, 'IntegratedTests.m'), ...
    'UseParallel', true);

disp(results);

if any([results.Failed])
    error( ...
        'GReS:IntegratedTestsFailed', ...
        'One or more integrated tests failed.');
end