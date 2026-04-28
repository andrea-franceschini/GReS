clc;
clear;
close all;

addpath('utils');

gresLog().setVerbosity(-2);

% Read input file
[meshProp, rock, blockSize, jointFamilies, druckerPrager, inputFiles] = ...
    readInputFile('config.xml');

% Import mesh data into the Mesh object
grid = Grid();
grid.importMesh(inputFiles.meshFile);

% Set parameters of the simulation
simParam = SimulationParameters(inputFiles.simParam);

% Create an object of the Materials class and read the materials file
mat = Materials(inputFiles.materials);

% Creating boundaries conditions.
bound = Boundaries(grid,inputFiles.boundaries);

% create the Discretizer (key-value pair input)
domain = Discretizer('grid',grid,...
                     'materials',mat,...
                     'boundaries',bound);

domain.addPhysicsSolver('Poromechanics');

% we keep in memory the initial state for repeated simulations
initState = copy(getState(domain));

% Micro to meso scale
rng(42);   % Set random seed
stiffnesses = createMat(domain.grid, meshProp, rock, blockSize, jointFamilies);

% Set materials
domain.materials.solid{1}.ConstLaw.setStiffnesses(stiffnesses);

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme. 
solver = NonLinearImplicit('simulationparameters', simParam, ...
                           'domains', domain);

fea.initState = initState;
fea.solver = solver;

incr = 9;
% incr = 3;
tol = 1e-2;
stepX = 1:incr:26;
nStep = length(stepX);

avgPs = zeros(nStep,1);
avgQs = zeros(nStep,1);

for id_forceX = 1:nStep
    forceX = stepX(id_forceX);
    [maxF, avgP, avgQ] = solvePQ(fea, forceX, druckerPrager, tol);
    fprintf('%d  %e  %e\n', forceX, avgP, avgQ);
    avgPs(id_forceX) = avgP;
    avgQs(id_forceX) = avgQ;
end

fid = fopen('Output/averages.txt','w');
for i = 1:length(avgPs)
    fprintf(fid,'%15.6e%15.6e\n', avgPs(i), avgQs(i));
end
fclose(fid);

[A, B] = fitCurve(avgPs, avgQs);

plotPQ(avgPs, avgQs, A, B, druckerPrager, 'Output/pq_plot.pdf');
