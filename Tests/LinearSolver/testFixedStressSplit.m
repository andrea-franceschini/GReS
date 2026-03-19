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

nn = [8,12,16];
linsolverTime = zeros(3,length(nn));
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

   F = -10; % vertical force
   state = applyMandelIC(domain.state,mat,mesh,F);

   % Solve the fully coupled Mandel Biot problem with matlab
   solver = NonLinearImplicit('simulationparameters',simParam,...
                              'domains',domain,...
                              'output',printUtils);
   solver.simulationLoop();

   % Get the time needed for the solve
   linsolverTime(1,i) = solver.linsolver.aTimeSolve;

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

   F = -10; % vertical force
   state = applyMandelIC(domain.state,mat,mesh,F);

   % Solve the Mandel Biot problem with matlab using the fixed stress split algorithm
   solver = FixedStressSplit('simulationparameters',simParam,...
                              'domains',domain,...
                              'maxiterations',30,...
                              'output',printUtils);
   solver.simulationLoop();

   % Sum the time contributions coming from the flow and mechanics solve steps
   linsolverTime(2,i) = solver.solverMech.aTimeSolve+solver.solverFlow.aTimeSolve;

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

   F = -10; % vertical force
   state = applyMandelIC(domain.state,mat,mesh,F);

   % Solve the Mandel Biot problem with the preconditioned iterative solver 
   % using the fixed stress split algorithm
   solver = FixedStressSplit('simulationparameters',simParam,...
                              'domains',domain,...
                              'maxiterations',30,...
                              'output',printUtils);
   solver.simulationLoop();

   % Sum the time contributions coming from the flow and mechanics solver 
   % (both solve and preconditioner computation steps)
   linsolverTime(3,i) = solver.solverMech.aTimeComp+solver.solverMech.aTimeSolve;
   linsolverTime(3,i) = linsolverTime(3,i) + solver.solverFlow.aTimeComp+solver.solverFlow.aTimeSolve;

end





%% Plot the speed comparison
figure;
plot(nn,linsolverTime(1,:),'ro-');
hold on
grid on
plot(nn,linsolverTime(2,:),'ko-');
plot(nn,linsolverTime(3,:),'go-');
hold off
legend('FC MATLAB','FS MATLAB','FS Iterative',Location='northwest')
xlabel('Number of Nodes per dimension');
ylabel('Computation Time (s)');
title('Solver Time Comparison');



