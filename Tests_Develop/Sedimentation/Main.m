close all;
% clear;
input_dir = 'Inputs/';
file_SimP = fullfile(input_dir,'simparam.xml');
file_Mat = fullfile(input_dir,'materials.xml');
file_Solver = fullfile(input_dir,'solver.xml');
% profile on;

%% ------------------------------ Set up the Domain -----------------------
% Set the simulation parameters for the non-linear solver.
simParam = SimulationParameters(file_SimP);

% Create an object of the Materials class and read the materials file
mat = Materials(file_Mat);

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Materials',mat);

domain.addPhysicsSolver(file_Solver);

% set initial conditions directly modifying the state object
phy = domain.physicsSolvers(domain.solverNames).grid;
z = phy.getCoordCenter(phy.getActiveDofs);
gamma_w = getFluid(mat).getFluidSpecWeight();
wLev = 1.; % level of the water table
domain.state.data.pressure = gamma_w*(wLev-z(:,3));

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme.
Solver = FCSolver(simParam,domain);

% Solve the problem
[simState] = Solver.NonLinearLoop();

domain.outstate.finalize()

% profile off
% profile viewer