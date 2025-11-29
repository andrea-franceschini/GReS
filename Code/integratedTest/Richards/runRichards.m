% Create and set the print utility
printUtils = OutState(topology,'output.xml');

% Create an object of the Materials class and read the materials file
mat = Materials(listMat(sim));

% Create object handling construction of Jacobian and rhs of the model
% linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);
domain = Discretizer('Grid',grid,...
                     'Materials',mat,...
                     'Boundaries',bound,...
                     'OutState',printUtils);

domain.addPhysicsSolver('solver.xml');

% set initial conditions directly modifying the state object
domain.state.data.pressure = getFluid(mat).getFluidSpecWeight()*(wLev-z);

% Solve the problem
Solver = FCSolver(simParam,domain);
[simState] = Solver.NonLinearLoop();

printUtils.finalize();