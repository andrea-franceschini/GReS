clear
close all

scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

cd(scriptDir);

fname = 'constantSlidingEFEM.xml';

simparams = SimulationParameters(fname);


% b = BlockStructuredMesh([0.0, 2.0; 0.0,0.5; 0.0, 4.0],[8,2,16],1);
% mesh = b.processGeometry();

mesh = structuredMesh(8,2,16,[0 2],[0, 0.5],[0 4]);


elems = Elements(mesh,2);
faces = Faces(mesh);
grid = struct('topology',mesh,'cells',elems,'faces',faces);
mat = Materials(fname);


% printUtils = OutState(mesh,"folderName",strcat("Output/ConstantSlidingEFEM"),"timeList",1,...
%                        "writeVtk",1,"flagMatFile",1,"matFileName",strcat("Output/ConstantSlidingHistory"));

printUtils = OutState(mesh);

bc = Boundaries(fname,grid);
% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Boundaries',bc,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);


domain.addPhysicsSolver(fname);


solver = GeneralSolver(simparams,domain);

solver.NonLinearLoop();
solver.finalizeOutput();

% get tangential gap
gt = abs(domain.state.data.fractureJump(2:3:end));
anGt = 0.1*sqrt(2);
tol = 1e-6;
assert(all(abs(gt - anGt)<tol),"Analytical solution is not matched")




% 

