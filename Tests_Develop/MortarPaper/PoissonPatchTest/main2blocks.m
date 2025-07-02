% patch test to check violation for RBF vs EB and SB 
% 2 simple blocks: 3x3 vs 4x4  
close all;
clear;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
cd(scriptDir);

% Set physical models 
model = ModelType("Poromechanics_FEM");

% Set parameters of the simulation
fileName = "simParam.dat";
simParam = SimulationParameters(fileName,model);

% Create the Mesh object
topology = Mesh();

% Set the mesh input file name
fileName = 'Mesh/cubeHexa8.vtk';
% Import the mesh data into the Mesh object
topology.importMesh(fileName);

% Create an object of the Materials class and read the materials file
fileName = 'materialsList.dat';
mat = Materials(model,fileName);


% Create an object of the "Elements" class and process the element properties
ngp = 2;
elems = Elements(topology,ngp);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);
%

% Degree of freedom manager 
%fname = 'dof.dat';
dofmanager = DoFManager(topology,model);

% Create object handling construction of Jacobian and rhs of the model
linSyst = Discretizer(model,simParam,dofmanager,grid,mat);

% Build a structure storing variable fields at each time step
linSyst.setState();

% Create and set the print utility
printUtils = OutState(model,topology,'outTime.dat','folderName','Output_PatchTest');


% % Write BC files programmatically with function utility 
% F = -10; % vertical force
% 
% % Collect BC input file in a list
% fileName = ["BCs/dirFlowTop.dat","BCs/neuPorotop.dat",...
%    "BCs/dirPoroLatY.dat","BCs/dirPoroLatX.dat","BCs/dirPoroBottom.dat"];
% %
writeBCfiles('BCs/fixBot','SurfBC','Dir',{'Poromechanics','x','y','z'},'bottom_fixed',0,0,topology,2); 
writeBCfiles('BCs/topLoad','SurfBC','Neu',{'Poromechanics','z'},'top_load',0,-1,topology,1); % left block lateral fix

fileName = ["BCs/fixBot.dat","BCs/topLoad.dat"];
% Create an object of the "Boundaries" class 
bound = Boundaries(fileName,model,grid);

% Print model initial state
printState(printUtils,linSyst);

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme. 
% Here, a built-in fully implict solution scheme is adopted with class
% FCSolver. This could be simply be replaced by a user defined function
Solver = FCSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,linSyst);
%
% Solve the problem
[simState] = Solver.NonLinearLoop();
%
% Finalize the print utility
printUtils.finalize()

