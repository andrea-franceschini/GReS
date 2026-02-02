clear
clc

tWall = tic();      % start wall-clock timer
tCPU  = cputime;    % start CPU timer

gresLog().setVerbosity(-1);

testFiles = {
    fullfile('Terzaghi','testTerzaghi.m')
    fullfile('SubDomains','testSubDomains.m')
    fullfile('MortarConvergence','testMortarPoisson.m')
    fullfile('Richards','testRichards.m')
    fullfile('ConstantSliding','testConstantSliding.m')
    fullfile('SingleCrackCompressed','testSingleCrackCompressed.m')
    fullfile('ConstantSlidingEFEM','testConstantSlidingEFEM.m')
};

results = runtests(testFiles);

elapsedWall = toc(tWall);
elapsedCPU  = cputime - tCPU;

fprintf("Elapsed wall-clock time: %1.2f s\n", elapsedWall);
fprintf("Elapsed CPU time:        %1.2f s\n", elapsedCPU);

if any([results.Failed])
    error("Some tests did not pass");
else
    disp("All tests ran successfully");
end
