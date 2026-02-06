% Create and set the print utility
printUtils = OutState('Input/output.xml');

% Create an object of the Materials class and read the materials file
mat = Materials(listMat(sim));

% Create object handling construction of Jacobian and rhs of the model
% linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);
domain = Discretizer('Grid',grid,...
                     'Materials',mat,...
                     'Boundaries',bound);

domain.addPhysicsSolver('Input/solver.xml');

% set initial conditions directly modifying the state object
domain.state.data.pressure = getFluid(mat).getSpecificWeight()*(wLev-z);

% Solve the problem
solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils);
solver.simulationLoop();

printUtils.finalize();