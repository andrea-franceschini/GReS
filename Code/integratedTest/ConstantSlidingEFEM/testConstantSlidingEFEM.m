clear
close all

scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

cd(scriptDir);

fname = 'constantSlidingEFEM.xml';

params = readInput(fname);

simparams = SimulationParameters(params.SimulationParameters);


mesh = structuredMesh(8,2,16,[0 2],[0, 0.5],[0 4]);


elems = Elements(mesh,2);
faces = Faces(mesh);
grid = struct('topology',mesh,'cells',elems,'faces',faces);
mat = Materials(params.Materials);


bc = Boundaries(grid,params.BoundaryConditions);
% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Boundaries',bc,...
  'Materials',mat,...
  'Grid',grid);


domain.addPhysicsSolvers(params.Solver);


solver = NonLinearImplicit('simulationparameters',simparams,...
  'domains',domain);
solver.simulationLoop();

% get tangential gap
gt = abs(domain.state.data.fractureJump(2:3:end));
anGt = 0.1*sqrt(2);
tol = 1e-6;
assert(all(abs(gt - anGt)<tol),"Analytical solution is not matched")




% 

