clear
close all

scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

cd(scriptDir);

input = readstruct('constantSliding.xml',AttributeSuffix="");

simparams = SimulationParameters(input.SimulationParameters);

[domains,interfaces] = buildModel('constantSliding.xml');

solver = NonLinearImplicit('simulationparameters',simparams,...
                           'domains',domains,...
                           'interface',interfaces);
solver.simulationLoop();

% get tangential gap
state = getState(interfaces{1});
gt = state.tangentialGap(1:2:end);
anGt = 0.1*sqrt(2);
tol = 1e-6;
assert(all(abs(gt - anGt)<tol),"Analytical solution is not matched")




% 

