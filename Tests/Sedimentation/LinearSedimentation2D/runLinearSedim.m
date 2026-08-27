close all;
% clear;
input_dir = 'Input/';
file_Mat = fullfile(input_dir,'materials.xml');
file_Solver = fullfile(input_dir,'solver.xml');

%% ------------------------------ Set up the Domain -----------------------
% Set the simulation parameters for the non-linear solver.
simParam = SimulationParameters('Start',0.,'End',1.0e3,...
      'DtInit',1e-0,'DtMin',1e-2,'DtMax',1e-0,'MaxNLIteration',30);

% Create an object of the Materials class and read the materials file
mat = Materials(file_Mat);
printUtils = OutState('outputFile','Outputs/Results2','printTimes',0:50:1000,"vtkFormat","ascii");

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Materials',mat);
domain.addPhysicsSolvers(file_Solver);

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme.
solver = EvolvingGrid('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils);

solver.simulationLoop();