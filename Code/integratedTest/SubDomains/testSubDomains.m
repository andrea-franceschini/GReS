close all;
clear;

profile on
scriptFullPath = mfilename('fullpath');
scriptDir = fileparts(scriptFullPath);
cd(scriptDir);
warning('off','MATLAB:nearlySingularMatrix');

% shortcut to define a model using a unique xml file
% useful when dealing with many domains
simparams = SimulationParameters('deepAquifer.xml');
output = OutState('deepAquifer.xml');
domain = buildModel('deepAquifer.xml');

% perform a fully coupled simulation
solver = NonLinearImplicit('simulationparameters',simparams,...
                           'domains',domain,...
                           'output',output);
solver.simulationLoop();

ref = load('referenceSol.mat');
comp = domain.outstate.matFile;
%%
tol = 1e-3;
assert(norm(ref.expDispl(:,1:end)-[comp.displacements])<tol)
assert(norm(ref.expPress(:,1:end)-[comp.pressure])<tol)


