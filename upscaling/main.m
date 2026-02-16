clc;
clear;
close all;

gresLog().setVerbosity(0);

% Read input file
[meshProp, rock, blockSize, jointFamilies, druckerPrager, inputFiles] = ...
    readInputFile('config.xml');

% Import mesh data into the Mesh object
mesh = Mesh();
mesh.importMesh(inputFiles.meshFile);
topology = Mesh();

% Set parameters of the simulation
simParam = SimulationParameters(inputFiles.simParam);

% Create an object of the Materials class and read the materials file
mat = Materials(inputFiles.materials);

% now the faces field can be omitted completely
grid = struct('topology',mesh,'cells',Elements(mesh));

% Creating boundaries conditions.
bound = Boundaries(inputFiles.boundaries,grid);

% create the Discretizer (key-value pair input)
domain = Discretizer('grid',grid,...
                     'materials',mat,...
                     'boundaries',bound);

domain.addPhysicsSolver(inputFiles.solver);

% we keep in memory the initial state for repeated simulations
initState = copy(getState(domain));

% Micro to meso scale
rng(42);   % Set random seed
stiffnesses = createMat(domain.grid, meshProp, rock, blockSize, jointFamilies);

% Set materials
domain.materials.db{1}.ConstLaw.setStiffnesses(stiffnesses);

%printUtils = OutState(fullfile('Input','output.xml'));

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme. 
solver = NonLinearImplicit('simulationparameters', simParam, ...
                           'domains', domain);%, ...
                           %'output', printUtils);

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

fid = fopen('averages.txt','w');
for i = 1:length(avgPs)
    fprintf(fid,'%15.6e%15.6e\n', avgPs(i), avgQs(i));
end
fclose(fid);

fitCurve(avgPs, avgQs);
