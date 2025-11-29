clear
close all

scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

cd(scriptDir);

simparams = SimulationParameters('constantSliding.xml');

[domains,interfaces] = buildModel('constantSliding.xml');

solver = ActiveSetContactSolver(simparams,domains,interfaces,5);

%solver.simParameters.setBackstepSkipFlag(1);

solver.NonLinearLoop();
solver.finalizeOutput();

% get tangential gap
gt = interfaces{1}.tangentialGap.curr;
anGt = 0.1*sqrt(2);
tol = 1e-6;
assert(all(abs(gt - anGt)<tol),"Analytical solution is not matched")




% 

