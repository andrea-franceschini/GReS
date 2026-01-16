clear
close all

scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

cd(scriptDir);

simparams = SimulationParameters('constantSliding.xml');

[domains,interfaces] = buildModel('constantSliding.xml');

solver = GeneralSolver(simparams,domains,interfaces);

solver.NonLinearLoop();
solver.finalizeOutput();

% get tangential gap
gt = interfaces{1}.state.tangentialGap(1:2:end);
anGt = 0.1*sqrt(2);
tol = 1e-6;
assert(all(abs(gt - anGt)<tol),"Analytical solution is not matched")




% 

