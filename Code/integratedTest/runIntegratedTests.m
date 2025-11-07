clear
clc


testFiles = {fullfile('Terzaghi','testTerzaghi.m');...
             fullfile('SubDomains','testSubDomains.m');...
             fullfile('MortarConvergence','testMortarPoisson.m');...
             fullfile('Richards','testRichards.m')
             fullfile('ConstantSliding','testConstantSliding.m')
             fullfile('SingleCrackCompressed','testSingleCrackCompressed.m')};

results = runtests(testFiles);

if any([results.Failed])
  error("Some test not passed");
else
  disp("All test run successfully")
end
