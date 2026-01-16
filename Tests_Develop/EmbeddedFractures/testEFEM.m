% simple 3D patch test to check minimal kernel in 3D

clear
close all
clc

%file = "simpleMechanics.xml";
file = "simpleMechanics.xml";

m = BlockStructuredMesh([0,1;0 1;0 1],[10,10,5],1);
mesh = processGeometry(m);

simParam = SimulationParameters(file);

elems = Elements(mesh,3);
faces = Faces(mesh);
grid = struct('topology',mesh,'cells',elems,'faces',faces);
mat = Materials(file);

% Create and set the print utility
printUtils = OutState(mesh,file);

% Create an object of the "Boundaries" class
bcs = Boundaries(file,grid);

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Boundaries',bcs,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);

domain.addPhysicsSolver(file);

solver = GeneralSolver(simParam,domain);
solver.NonLinearLoop();
domain.outstate.finalize();

% finalize EFEM outstate
%getPhysicsSolver(domain,"EmbeddedFractureMechanics").outFracture.finalize();

