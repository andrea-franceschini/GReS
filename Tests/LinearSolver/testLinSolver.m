% This script runs a series of simulation benchmarks comparing two solver
% configurations (MATLAB direct and iterative) for flow and poromechanics
% problems. For each problem type it:
%   - reads simulation and material input files,
%   - builds structured meshes of increasing resolution,
%   - constructs discretization, boundary and material objects,
%   - runs the NonLinearImplicit solver,
%   - records linear solver computation times,
% and finally plots solver time comparisons for flow and mechanics tests.

close all;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
cd(scriptDir);

gresLog().setVerbosity(2)

%% SINGLE PHYSICS - SINGLE PHASE FLOW MATLAB
clear
clc

fileName = "flow/flowMatlab.xml";
params = readInput(fileName);

% Set parameters of the simulation
simParam = SimulationParameters(params.SimulationParameters);

% Create an object of the Materials class and read the materials file
mat = Materials(params.Materials);


% Number of nodes for one dimension
nn = [16,25,35,45,55];
linsolverTime = zeros(2,length(nn));
for i = 1:length(nn)
   mesh = structuredMesh(nn(i),nn(i),nn(i),[0 100],[0 100],[0 10]);

   gaussOrder = 2;
   elems = Elements(mesh,gaussOrder);
   faces = Faces(mesh);
   gridd = struct('topology',mesh,'cells',elems,'faces',faces);


   bound = Boundaries(gridd,params.BoundaryConditions);

   printUtils = OutState(params.Output);

   domain = Discretizer('Boundaries',bound,...
                        'Materials',mat,...
                        'Grid',gridd);
   domain.addPhysicsSolvers(params.Solver);

   solver = NonLinearImplicit('simulationparameters',simParam,...
                              'domains',domain,...
                              'output',printUtils);
   solver.simulationLoop();

   % Get the time needed for the solve step
   linsolverTime(1,i) = solver.linsolver.aTimeSolve;

end






%% SINGLE PHYSICS - SINGLE PHASE FLOW ITERATIVE SOLVER
clc

fileName = "flow/flowIterative.xml";
params = readInput(fileName);

% Set parameters of the simulation
simParam = SimulationParameters(params.SimulationParameters);

% Create an object of the Materials class and read the materials file
mat = Materials(params.Materials);

for i = 1:length(nn)
   mesh = structuredMesh(nn(i),nn(i),nn(i),[0 100],[0 100],[0 10]);

   gaussOrder = 2;
   elems = Elements(mesh,gaussOrder);
   faces = Faces(mesh);
   gridd = struct('topology',mesh,'cells',elems,'faces',faces);


   bound = Boundaries(gridd,params.BoundaryConditions);

   printUtils = OutState(params.Output);

   domain = Discretizer('Boundaries',bound,...
                        'Materials',mat,...
                        'Grid',gridd);
   domain.addPhysicsSolvers(params.Solver);

   solver = NonLinearImplicit('simulationparameters',simParam,...
                              'domains',domain,...
                              'output',printUtils);
   solver.simulationLoop();

   % Get the time needed for the solve step and the preconditioner computation time
   linsolverTime(2,i) = solver.linsolver.aTimeComp+solver.linsolver.aTimeSolve;

end

%%
% Plot the speed comparison
figure;
plot(nn,linsolverTime(1,:),'ro-');
hold on
grid on
plot(nn,linsolverTime(2,:),'ko-');
hold off
legend("MATLAB","Iterative",Location="northwest")
xlabel('Number of Nodes per dimension');
ylabel('Computation Time (s)');
title('Solver Time Comparison Flow');



%% SINGLE PHYSICS - POROMECHANICS MATLAB
clear
clc

fileName = "mechanics/mechMatlab.xml";
params = readInput(fileName);

% Set parameters of the simulation
simParam = SimulationParameters(params.SimulationParameters);

% Create an object of the Materials class and read the materials file
mat = Materials(params.Materials);


% Number of nodes for one dimension
nn = [8,12,16,20,24];
linsolverTime = zeros(2,length(nn));
for i = 1:length(nn)
   mesh = structuredMesh(nn(i),nn(i),nn(i),[0 100],[0 100],[0 10]);

   gaussOrder = 2;
   elems = Elements(mesh,gaussOrder);
   faces = Faces(mesh);
   gridd = struct('topology',mesh,'cells',elems,'faces',faces);


   bound = Boundaries(gridd,params.BoundaryConditions);

   printUtils = OutState(params.Output);

   domain = Discretizer('Boundaries',bound,...
                        'Materials',mat,...
                        'Grid',gridd);
   domain.addPhysicsSolvers(params.Solver);

   solver = NonLinearImplicit('simulationparameters',simParam,...
                              'domains',domain,...
                              'output',printUtils);
   solver.simulationLoop();

   % Get the time needed for the solve step 
   linsolverTime(1,i) = solver.linsolver.aTimeSolve;

end






%% SINGLE PHYSICS - POROMECHANICS ITERATIVE SOLVER
clc

fileName = "mechanics/mechIterative.xml";
params = readInput(fileName);

% Set parameters of the simulation
simParam = SimulationParameters(params.SimulationParameters);

% Create an object of the Materials class and read the materials file
mat = Materials(params.Materials);

for i = 1:length(nn)
   mesh = structuredMesh(nn(i),nn(i),nn(i),[0 100],[0 100],[0 10]);

   gaussOrder = 2;
   elems = Elements(mesh,gaussOrder);
   faces = Faces(mesh);
   gridd = struct('topology',mesh,'cells',elems,'faces',faces);


   bound = Boundaries(gridd,params.BoundaryConditions);

   printUtils = OutState(params.Output);

   domain = Discretizer('Boundaries',bound,...
                        'Materials',mat,...
                        'Grid',gridd);
   domain.addPhysicsSolvers(params.Solver);

   solver = NonLinearImplicit('simulationparameters',simParam,...
                              'domains',domain,...
                              'output',printUtils);
   solver.simulationLoop();

   % Get the time needed for the solve step and the preconditioner computation time
   linsolverTime(2,i) = solver.linsolver.aTimeComp+solver.linsolver.aTimeSolve;

end

%%
% Plot the speed comparison
figure
plot(nn,linsolverTime(1,:),'ro-');
hold on
grid on
plot(nn,linsolverTime(2,:),'ko-');
hold off
legend("MATLAB","Iterative",Location="northwest")
xlabel('Number of Nodes per dimension');
ylabel('Computation Time (s)');
title('Solver Time Comparison Mechanics');
