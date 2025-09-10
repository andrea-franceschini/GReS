close all;
% clear;
input_dir = 'Inputs';
output_dir = 'Outputs';
figures_dir = 'Figs';

%% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("VariabSatFlow_FVTPFA");

%% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = fullfile(input_dir,'simParam.dat');
simParam = SimulationParameters(fileName,model);

%% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();

% Set the input file name
% fileName = fullfile(input_dir,'Mesh','Column30.msh');
fileName = fullfile(input_dir,'Mesh','background.vtk');

% Import mesh data into the Mesh object
% topology.importGMSHmesh(fileName);
topology.importMesh(fileName);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model,topology);

% spot cells that do not have any neighboring cells (no internal faces)
intFaces = all(faces.faceNeighbors,2);
validCells = unique(faces.faceNeighbors(intFaces,:));
id = true(topology.nCells,1);
id(validCells) = false;

% remove invalid cells from fluidMesh
topology.cells(id,:) = [];
topology.cellTag(id) = [];
topology.nCells = topology.nCells - sum(id);
topology.cellVTKType(id) = [];

%reprocess faces
faces = Faces(model,topology);
intFaces = all(faces.faceNeighbors,2);
validCells = unique(faces.faceNeighbors(intFaces,:));
id = true(topology.nCells,1);
id(validCells) = false;
assert(sum(id)==0,'something wrong')

topology = setSurfacesFluid(topology,faces);

% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,2);

% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);

%% ----------------------------- MATERIALS -----------------------------
% Set the input file name
fileName = fullfile(input_dir,'materialsList.dat');

% Create an object of the Materials class and read the materials file
mat = Materials(model,fileName);

%% ----------------------- DOF Manager -----------------------------
% Degree of freedom manager 
%fname = 'dof.dat';
dofmanager = DoFManager(topology,model);

% Create and set the print utility
printUtils = OutState(model,topology,fullfile(input_dir,'outTime.dat'), ...
    'folderName','Outputs','flagMatFile',true);

%% ----------------------- Boundary Condition -----------------------------
% Creating and Appling boundaries conditions.
% fileName = ["Inputs/BC_Bottom.dat","Inputs/BC_Top.dat"];
fileName = ["BCs/press_top.dat","BCs/press_bot.dat"];
bound = Boundaries(fileName,model,grid);

%% ----------------------- Discretizer -----------------------------
% Create object handling construction of Jacobian and rhs of the model
% linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);
domain = Discretizer('ModelType',model,...
                     'SimulationParameters',simParam,...
                     'DoFManager',dofmanager,...
                     'Boundaries',bound,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);

%% ----------------------- Initial Condition -----------------------------
% Build a structure storing variable fields at each time step
% state = linSyst.setState();

% set initial conditions directly modifying the state object
domain.state.data.pressure(:) = -10*9.8066e3;

%% ----------------------- Solver -----------------------------
% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme. 
% Here, a built-in fully implict solution scheme is adopted with class
% FCSolver. This could be simply be replaced by a user defined function
Solver = FCSolver(domain,'SaveRelError',true,'SaveBStepInf',true);

% Solve the problem
[simState] = Solver.NonLinearLoop();

% Finalize the print utility
printUtils.finalize()










%% ----------------------- Extra -----------------------------

function msh = setSurfacesFluid(msh,faces)
% get external faces
id = find(~all(faces.faceNeighbors,2));

ids = [faces.mapN2F(id) faces.mapN2F(id)+1 faces.mapN2F(id)+2 faces.mapN2F(id)+3];
msh.surfaces = faces.nodes2Faces(ids);

nS = size(msh.surfaces,1);
% get surface on top
msh.surfaceTag = ones(nS,1);
msh.nSurfaces = nS;

% find surfaces on top
tol = 1e-4;
cZ = msh.coordinates(msh.surfaces',3);
cZ = reshape(cZ,4,[]);
minZ = min(msh.coordinates(:,3),[],'all');
maxZ = max(msh.coordinates(:,3),[],'all');
idTop = all(abs(cZ-maxZ)<tol,1);
idBot = all(abs(cZ-minZ)<tol,1);
msh.surfaceTag(idBot) = 2;
msh.surfaceTag(idTop) = 3;
msh.nSurfaceTag = 3;
msh.surfaceVTKType = 9*ones(nS,1);
msh.surfaceNumVerts = 4*ones(nS,1);


end