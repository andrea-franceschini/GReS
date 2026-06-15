% Standalone driver script for a vertically compressed elastoplastic cube
clear; clc;

simParam = SimulationParameters(fullfile('Input','simparam.xml'));

% 1. Create a structured 3D Cube Mesh (1m x 1m x 1m)
X = 5; Y = 1; Z = 10;
nx = 10; ny = 5; nz = 10;
grid = structuredMesh(nx,ny,nz,[0 X],[0 Y],[0 Z]);

material = Materials(fullfile('Input','materials.xml'));

bound = Boundaries(grid,fullfile('Input','boundaries.xml'));

printUtils = OutState('printTimes',0:1:20,'outputFile',"Output/results");

domain = Discretizer('Grid',grid,...
                     'Materials',material,...
                     'Boundaries',bound);

domain.addPhysicsSolver('Poromechanics');

solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domain, ...
                           'output',printUtils);
solver.simulationLoop();