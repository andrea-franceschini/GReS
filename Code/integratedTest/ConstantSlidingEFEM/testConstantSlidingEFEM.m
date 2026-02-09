clear
close all

scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

cd(scriptDir);

fname = 'constantSlidingEFEM.xml';

simparams = SimulationParameters(fname);


mesh = structuredMesh(8,2,16,[0 2],[0, 0.5],[0 4]);


elems = Elements(mesh,2);
faces = Faces(mesh);
grid = struct('topology',mesh,'cells',elems,'faces',faces);
mat = Materials(fname);


bc = Boundaries(fname,grid);
% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Boundaries',bc,...
                     'Materials',mat,...
                     'Grid',grid);


domain.addPhysicsSolver(fname);


solver = NonLinearImplicit('simulationparameters',simparams,...
                           'domains',domain);
solver.simulationLoop();

% get tangential gap
gt = abs(domain.state.data.fractureJump(2:3:end));
anGt = 0.1*sqrt(2);
tol = 1e-6;
assert(all(abs(gt - anGt)<tol),"Analytical solution is not matched")




% 

