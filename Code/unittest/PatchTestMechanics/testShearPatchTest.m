clear
close all

scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

cd(scriptDir);

fname = 'shearPatch.xml';

simparams = SimulationParameters(fname);


b = BlockStructuredMesh([0.0, 1.0;0.0 1.0; 0.0 1.0],[1,1,1],1);
mesh = b.processGeometry();


elems = Elements(mesh,2);
faces = Faces(mesh);
grid = struct('topology',mesh,'cells',elems,'faces',faces);
mat = Materials(fname);


printUtils = OutState(mesh,"writeVtk",0,"flagMatFile",0);


bc = Boundaries(fname,grid);
% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Boundaries',bc,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);


domain.addPhysicsSolver(fname);


solver = GeneralSolver(simparams,domain);

solver.NonLinearLoop();




% 

