% RUNSOFTSOILCREEPVALIDATION - Model validation against laboratory triaxial/creep tests
% based on Isotton et al. (2019).
clear; clc;

% Load simulation parameters
simParam = SimulationParameters(fullfile('Input', 'simparam.xml'));

% 1. Create a representative 3D single-element mesh (1m x 1m x 1m)
X = 1.0; Y = 1.0; Z = 1.0;
nx = 10; ny = 10; nz = 10;
grid = structuredMesh(nx, ny, nz, [0 X], [0 Y], [0 Z]);

% 2. Load Soft Soil Creep Material Database
material = Materials(fullfile('Input', 'materials.xml'));

% 3. Set Boundary Conditions for Triaxial Compression / Creep Test
bound = Boundaries(grid, fullfile('Input', 'boundaries.xml'));

% 4. Output utility
printUtils = OutState('printTimes', 0:0.1:20.0, 'outputFile', "Output/resultsSoftSoilCreep", 'matFileName',"Output/results");

% 5. Assemble Discretizer Domain
domain = Discretizer('Grid', grid, ...
                     'Materials', material, ...
                     'Boundaries', bound);

domain.addPhysicsSolver('Poromechanics');

% Initialize initial stress state (compression-negative in GReS)
initialStress(:, 1:3) = -repmat([0.43 0.43 1.00],size(domain.getState.stress, 1),1);
initialStress(:, 4:6) = 0.0; % Zero shear stresses
domain.setState(initialStress, "stress");

% 6. Launch Non-Linear Implicit Simulation Loop
solver = NonLinearImplicit('simulationparameters', simParam, ...
                           'domains', domain, ...
                           'output', printUtils);
solver.simulationLoop();

fprintf('Soft Soil Creep validation simulation completed successfully.\n');

%%

% res = load("Output/results.mat");
% 
% out = res.output;
% % Extract results from the simulation
% stressResults = res.output(:).stress;
% % strainResults = res.strain(:,3);
% timeResults = res.output(:).time;