close all;

% Resolve simulation input files relative to this script directory
scriptDir = fileparts(mfilename('fullpath'));
inputDir = fullfile(scriptDir, 'Sedimentation2D', 'Input');
file_Mat = fullfile(inputDir, 'materials.xml');
file_Solver = fullfile(inputDir, 'solver.xml');

% Configure linear solver to force iterative solver
linParams = struct('useIterative', 1);

% Configure simulation time parameters with a drastically shortened horizon.
% With sedimentation rate 0.003 m/s and cell height threshold 0.10 m
simParam = SimulationParameters('Start', 0.0, 'End', 300.0, ...
   'DtInit', 10.0, 'DtMin', 1e-2, 'DtMax', 10.0, 'MaxNLIteration', 8, ...
   'LinearSolver', linParams);

% Instantiate materials database from XML configuration
mat = Materials(file_Mat);

% Read solver setup and inject absolute path to sedimentation history map
solverInput = readInput(file_Solver);
solverInput.Sedimentation.SedimentMap.file = fullfile(inputDir, 'mapsSedimentation.xml');

% Build spatial discretizer and register sedimentation physics solver
domain = Discretizer('Materials', mat);
domain.addPhysicsSolvers(solverInput);

% Construct evolving grid solution scheme with dummy output utility
printUtils = OutState();
solver = EvolvingGrid('simulationparameters', simParam, ...
                      'domains', domain, ...
                      'output', printUtils, ...
                      'freezeCbAt', 1e-1);

% Run nonlinear time-stepping loop
solver.simulationLoop();
solverSedim = solver;
