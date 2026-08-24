% This script runs a set of poromechanics simulations using three solver
% variants (fully coupled, fixed-stress split, and fixed-stress iterative)
% for increasing mesh resolutions and compares linear solver times.
%
% High-level steps:
% 1. Set working directory to the script location so relative file I/O works.
% 2. Configure logging verbosity.
% 3. For three different configurations (fully coupled, fixed-stress split,
%    and fixed-stress iterative) the script:
%    - Reads an XML input file with simulation parameters.
%    - Constructs material, mesh, element, face, boundary and discretization
%      objects required by the simulation framework.
%    - Applies Mandel-type initial conditions (applyMandelIC) including a
%      vertical force F.
%    - Creates the appropriate solver object (NonLinearImplicit or
%      FixedStressSplit) with the provided simulation parameters and runs
%      the time/iteration loop via simulationLoop().
%    - Records the linear solver assembly and solve times for later
%      comparison.
% 4. After running all simulations for the chosen mesh sizes, the script
%    plots a comparison of the accumulated linear solver times versus the
%    number of nodes per dimension for each solver type.
%
% Notes on variables:
% - nn: array of numbers of nodes per dimension used to build structured
%   meshes.
% - linsolverTime: 3-by-length(nn) matrix storing solver assembly+solve
%   times for each solver type and mesh resolution. Row 1 = fully coupled,
%   Row 2 = fixed-stress split, Row 3 = fixed-stress iterative.
% - The script relies on framework classes/functions: readInput,
%   SimulationParameters, Materials, structuredMesh, Elements, Faces,
%   Boundaries, OutState, Discretizer, applyMandelIC, NonLinearImplicit,
%   FixedStressSplit, and their documented properties (e.g. aTimeComp,
%   aTimeSolve).
%
% This comment block describes the intent and structure of the script and
% what each major section accomplishes.


close all;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
cd(scriptDir);

gresLog().setVerbosity(2)

%% COUPLED PHYSICS - FULLY COUPLED MATLAB
clc

fileName = 'fixedStress/fullyCoupledMatlab.xml';
params = readInput(fileName);

% Set parameters of the simulation
simParam = SimulationParameters(params.SimulationParameters);

% Create an object of the Materials class and read the materials file
mat = Materials(params.Materials);

nn = [8,12,16,20];
linsolverTime = zeros(4,length(nn));
totSimTime = zeros(4,length(nn));
for i = 1:length(nn)
   gridd = structuredMesh(nn(i),nn(i),nn(i),[0 100],[0 100],[0 10]);

   bound = Boundaries(gridd,params.BoundaryConditions);

   printUtils = OutState(params.Output);

   domain = Discretizer('Boundaries',bound,...
                        'Materials',mat,...
                        'Grid',gridd);
   domain.addPhysicsSolvers(params.Solver);

   F = -10; % vertical force
   state = applyMandelIC(domain.state,mat,gridd,F);

   % Solve the fully coupled Mandel Biot problem with matlab
   solver = NonLinearImplicit('simulationparameters',simParam,...
                              'domains',domain,...
                              'output',printUtils);
   t0 = tic;
   solver.simulationLoop();
   totSimTime(1,i) = toc(t0);

   % Get the time needed for the solve
   linsolverTime(1,i) = solver.linsolver.getTotalTime();

end

%% COUPLED PHYSICS - FULLY COUPLED ITERATIVE
clc

fileName = 'fixedStress/fullyCoupledIterative.xml';
params = readInput(fileName);

% Set parameters of the simulation
simParam = SimulationParameters(params.SimulationParameters);

% Create an object of the Materials class and read the materials file
mat = Materials(params.Materials);

for i = 1:length(nn)
   gridd = structuredMesh(nn(i),nn(i),nn(i),[0 100],[0 100],[0 10]);

   bound = Boundaries(gridd,params.BoundaryConditions);

   printUtils = OutState(params.Output);

   domain = Discretizer('Boundaries',bound,...
                        'Materials',mat,...
                        'Grid',gridd);
   domain.addPhysicsSolvers(params.Solver);

   F = -10; % vertical force
   state = applyMandelIC(domain.state,mat,gridd,F);

   % Solve the fully coupled Mandel Biot problem with matlab
   solver = NonLinearImplicit('simulationparameters',simParam,...
                              'domains',domain,...
                              'output',printUtils);
   t0 = tic;
   solver.simulationLoop();
   totSimTime(2,i) = toc(t0);

   % Get the time needed for the solve
   linsolverTime(2,i) = solver.linsolver.getTotalTime();
   
end





%% COUPLED PHYSICS - FIXED STRESS SPLIT MATLAB
clc

fileName = 'fixedStress/fixedStressMatlab.xml';
params = readInput(fileName);

% Set parameters of the simulation
simParam = SimulationParameters(params.SimulationParameters);

% Create an object of the Materials class and read the materials file
mat = Materials(params.Materials);

% Number of nodes for one dimension
for i = 1:length(nn)
   gridd = structuredMesh(nn(i),nn(i),nn(i),[0 100],[0 100],[0 10]);


   bound = Boundaries(gridd,params.BoundaryConditions);

   printUtils = OutState(params.Output);

   domain = Discretizer('Boundaries',bound,...
                        'Materials',mat,...
                        'Grid',gridd);
   domain.addPhysicsSolvers(params.Solver);

   F = -10; % vertical force
   state = applyMandelIC(domain.state,mat,gridd,F);

   % Solve the Mandel Biot problem with matlab using the fixed stress split algorithm
   solver = FixedStressSplit('simulationparameters',simParam,...
                              'domains',domain,...
                              'maxiterations',30,...
                              'output',printUtils);
   t0 = tic;
   solver.simulationLoop();
   totSimTime(3,i) = toc(t0);

   % Sum the time contributions coming from the flow and mechanics solve steps
   linsolverTime(3,i) = solver.solverMech.getTotalTime()+solver.solverFlow.getTotalTime();

end






%% COUPLED PHYSICS - FIXED STRESS SPLIT ITERATIVE
clc

fileName = 'fixedStress/fixedStressIterative.xml';
params = readInput(fileName);

% Set parameters of the simulation
simParam = SimulationParameters(params.SimulationParameters);

% Create an object of the Materials class and read the materials file
mat = Materials(params.Materials);

for i = 1:length(nn)
   gridd = structuredMesh(nn(i),nn(i),nn(i),[0 100],[0 100],[0 10]);


   bound = Boundaries(gridd,params.BoundaryConditions);

   printUtils = OutState(params.Output);

   domain = Discretizer('Boundaries',bound,...
                        'Materials',mat,...
                        'Grid',gridd);
   domain.addPhysicsSolvers(params.Solver);

   F = -10; % vertical force
   state = applyMandelIC(domain.state,mat,gridd,F);

   % Solve the Mandel Biot problem with the preconditioned iterative solver 
   % using the fixed stress split algorithm
   solver = FixedStressSplit('simulationparameters',simParam,...
                              'domains',domain,...
                              'maxiterations',30,...
                              'output',printUtils);
   t0 = tic;
   solver.simulationLoop();
   totSimTime(4,i) = toc(t0);

   % Sum the time contributions coming from the flow and mechanics solver 
   linsolverTime(4,i) = solver.solverMech.getTotalTime() + solver.solverFlow.getTotalTime();
   
end



%% Plot the speed comparison
figure;
plot(nn,linsolverTime(1,:),'ro-');
hold on
grid on
plot(nn,linsolverTime(2,:),'ko-');
plot(nn,linsolverTime(3,:),'go-');
plot(nn,linsolverTime(4,:),'bo-');
hold off
legend('FC MATLAB','FC Iterative','FS MATLAB','FS Iterative',Location='northwest')
xlabel('Number of Nodes per dimension');
ylabel('Computation Time (s)');
title('Solver Time Comparison');


figure;
plot(nn,totSimTime(1,:),'ro-');
hold on
grid on
plot(nn,totSimTime(2,:),'ko-');
plot(nn,totSimTime(3,:),'go-');
plot(nn,totSimTime(4,:),'bo-');
hold off
legend('FC MATLAB','FC Iterative','FS MATLAB','FS Iterative',Location='northwest')
xlabel('Number of Nodes per dimension');
ylabel('Simulation Time (s)');
title('Simulation Time Comparison');