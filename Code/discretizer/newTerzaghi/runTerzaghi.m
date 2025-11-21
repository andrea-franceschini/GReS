% model = ModelManager();
% model.createModel('terzaghi.xml');
% 
% 
% model.runProblem();

close all
clc

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

cd(scriptDir);


% new model without requiring ModelType and DoFManager. Input file for physicsSolver
% is required

simparams = SimulationParameters('simparam.xml');

mesh = Mesh();
mesh.importMesh('Mesh/Column_hexa.msh');

elems = Elements(mesh,2);
faces = Faces(mesh);

grid = struct('topology',mesh,'cells',elems,'faces',faces);

mat = Materials('materials.xml');

bc = Boundaries('boundaryConditions.xml',grid);

printUtils = OutState(mesh,'output.xml');


% create the Discretizer (key-value pair)
domain = DiscretizerNew('grid',grid,...
                        'materials',mat,...
                        'boundaries',bc,...
                        'outstate',printUtils);

domain.addPhysicsSolver('solver.xml');

% manually apply initial conditions
state = domain.getState();
applyTerzaghiIC(state,mat,mesh,-10);

solv = FCSolver(simparams,domain);

stat = solv.NonLinearLoop();

domain.outstate.finalize();
