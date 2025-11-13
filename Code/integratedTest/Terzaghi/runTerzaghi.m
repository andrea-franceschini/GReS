% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);
% Change the current directory to the script's directory
cd(scriptDir);
% Set parameters of the simulation
fileName = "simParam.xml";
simParam = SimulationParameters(fileName,model);
simParam.setVerbosity(0);
% Create the Mesh object
topology = Mesh();
% Set the mesh input file name
fileName = 'Mesh/Column_hexa.msh';
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);
% Create an object of the Materials class and read the materials file
fileName = 'materials.xml';
mat = Materials(fileName);
% Create an object of the "Elements" class and process the element properties
gaussOrder = 2;
elems = Elements(topology,gaussOrder);
% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);
%
grid = struct('topology',topology,'cells',elems,'faces',faces);
%
dofmanager = DoFManager(topology,model);
% Create and set the print utility
%printUtils = OutState(model,topology,'flagMatFile',true,'writeVtk',true,'timeList',[1,5,10]);
printUtils = OutState(model,topology,'output.xml');
%
F = -10; % vertical force
setTerzaghiBC('BCs',F,topology);

% Collect BC input file in a list
fileName = ["BCs/dirFlowTop.dat","BCs/neuPorotop.dat",...
   "BCs/dirPoroLatY.dat","BCs/dirPoroLatX.dat","BCs/dirPoroBottom.dat"];
%
% Create an object of the "Boundaries" class 
bound = Boundaries('boundaryConditions.xml',model,grid);

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('ModelType',model,...
                     'SimulationParameters',simParam,...
                     'DoFManager',dofmanager,...
                     'Boundaries',bound,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);
applyTerzaghiIC(domain.state,mat,topology,F);
printState(domain);
Solver = FCSolver(domain);

%calling analytical solution script
%Terzaghi_analytical(topology, mat, 10)

[simState] = Solver.NonLinearLoop();