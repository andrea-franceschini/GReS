
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

nn = [16];
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
   solver = NonLinearImplicit('simulationparameters',simParam,...
                              'domains',domain,...
                              'output',printUtils);
   solver.simulationLoop();
   linsolverTime(1,i) = solver.linsolver.aTimeComp+solver.linsolver.aTimeSolve;

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
   solver = FixedStressSplit('simulationparameters',simParam,...
                              'domains',domain,...
                              'maxiterations',30,...
                              'output',printUtils);
   solver.simulationLoop();
   linsolverTime(2,i) = solver.solverMech.aTimeComp+solver.solverMech.aTimeSolve;
   linsolverTime(2,i) = linsolverTime(2,i) + solver.solverFlow.aTimeComp+solver.solverFlow.aTimeSolve;

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
   solver = FixedStressSplit('simulationparameters',simParam,...
                              'domains',domain,...
                              'maxiterations',30,...
                              'output',printUtils);
   solver.simulationLoop();
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



