close all;
clear;

profile on
scriptFullPath = mfilename('fullpath');
scriptDir = fileparts(scriptFullPath);
cd(scriptDir);
warning('off','MATLAB:nearlySingularMatrix');

% shortcut to define a model using a unique xml file
% useful when dealing with many domains
domain = buildModel('domain.xml');

domain.simparams.setVerbosity(0);

% perform a fully coupled simulation
solver = FCSolver(domain);
[simState] = solver.NonLinearLoop();

ref = load('referenceSol.mat');
comp = domain.outstate.results;

tol = 1e-3;
assert(norm(ref.expDispl(:,2:end)-[comp.expDispl])<tol)
assert(norm(ref.expPress(:,2:end)-[comp.expPress])<tol)


