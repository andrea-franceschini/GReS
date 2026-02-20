% close all;
% clear;
input_dir = 'Input/';
file_SimP = fullfile(input_dir,'simparam.xml');
% file_Mat = fullfile(input_dir,'materialsElastic.xml');
 file_Mat = fullfile(input_dir,'materials.xml');
file_Output = fullfile(input_dir,'output.xml');
file_Solver = fullfile(input_dir,'solver.xml');

%% ------------------------------ Set up the Domain -----------------------
% Set the simulation parameters for the non-linear solver.
simParam = SimulationParameters(file_SimP);

% Create an object of the Materials class and read the materials file
mat = Materials(file_Mat);
printUtils = OutState(file_Output);

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Materials',mat);
domain.addPhysicsSolver(file_Solver);

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme.
solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils);
solver.simulationLoop();

% % % % %fprintf("Bottom pressure -"domain.state.data.pressure(1))
% % % % sed = domain.getPhysicsSolver("Sedimentation");
% % % % 
% % % % load('/scratches/flash/castro/GReS/Code/integratedTest/PillarSed/Outputs/Results.mat')
% % % % pos = 9;
% % % % comp=output(pos).compaction(1:4:end);
% % % % xpos = sed.grid.coordZ(1:length(comp));
% % % % plot(xpos,comp)
% % % % 
% % % % % def = domain.state.data.cellDefm(1:4:sed.grid.ndofs)
% % % % % pres = domain.state.data.pressure(1:4:sed.grid.ndofs)
% % % % % plot(xpos,abs(def))
% % % % hold on