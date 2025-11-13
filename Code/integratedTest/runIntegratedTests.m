clear
clc

tic
t0 = cputime;
testFiles = {fullfile('Terzaghi','testTerzaghi.m');...
             fullfile('SubDomains','testSubDomains.m');...
             fullfile('MortarConvergence','testMortarPoisson.m');...
             fullfile('Richards','testRichards.m')
             fullfile('ConstantSliding','testConstantSliding.m')
             fullfile('SingleCrackCompressed','testSingleCrackCompressed.m')
             };

results = runtests(testFiles);

t = toc;
fprintf("Elasped wall-clock time: %1.2f s \n", t)
fprintf("Elasped CPU time: %1.2f s \n", cputime-t0)

if any([results.Failed])
  error("Some test not passed");
else
  disp("All test run successfully")
end
