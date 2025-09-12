%% RUN SPFLOW PROBLEM ON THE STRUCTURED MESH

% define permeability field using cubic law d^3/12
fluidMesh = Mesh();
fluidMesh.importMesh('Mesh/background.vtk');

% get surfaces
% 1 - internal grain surfaces
% 2 - bottom surfaces
% 3 - top surfaces
%fluidMesh = setSurfacesFluid(fluidMesh);

% run flow model
model = ModelType("VariabSatFlow_FVTPFA");
simParam = SimulationParameters('simParam.xml');
fprintf('Processing faces \n')
faces = Faces(model,fluidMesh);

% spot cells that do not have any neighboring cells (no internal faces)
intFaces = all(faces.faceNeighbors,2);
validCells = unique(faces.faceNeighbors(intFaces,:));
id = true(fluidMesh.nCells,1);
id(validCells) = false;

% remove invalid cells from fluidMesh
fluidMesh.cells(id,:) = [];
fluidMesh.cellTag(id) = [];
fluidMesh.nCells = fluidMesh.nCells - sum(id);
fluidMesh.cellVTKType(id) = [];

%reprocess faces
faces = Faces(model,fluidMesh);
intFaces = all(faces.faceNeighbors,2);
validCells = unique(faces.faceNeighbors(intFaces,:));
id = true(fluidMesh.nCells,1);
id(validCells) = false;
assert(sum(id)==0,'something wrong')

fluidMesh = setSurfacesFluid(fluidMesh,faces);
fprintf('Processing elements \n')
elems = Elements(fluidMesh,2);

gridFluid = struct('topology',fluidMesh,'cells',elems,'faces',faces);
material = Materials(model,'Materials/materialsVSFlow.dat');
dofFluid = DoFManager(fluidMesh,model);

writeBCfiles('BCs/press_top','SurfBC','Dir',"VariablySaturatedFlow",'Top_pressure',0,-7.354950e+03,fluidMesh,3);
writeBCfiles('BCs/press_bot','SurfBC','Dir',"VariablySaturatedFlow",'Bottom_pressure',0,-98066,fluidMesh,2);

bc = Boundaries(["BCs/press_top.dat","BCs/press_bot.dat"],model,gridFluid);
printUtils = OutState(model,fluidMesh,'outTime.dat','writeVtk',true,'folderName','Output_flow');
linSystFluid = Discretizer('ModelType',model,...
                     'SimulationParameters',simParam,...
                     'DoFManager',dofFluid,...
                     'Boundaries',bc,...
                     'OutState',printUtils,...
                     'Materials',material,...
                     'Grid',gridFluid);
linSystFluid.state.data.pressure(:) = -10*9.8066e3;

% solverFlow = FCSolver(model,simParam,dofFluid,gridFluid,material,bc,printUtils,linSystFluid);
solverFlow = FCSolver(linSystFluid,'SaveRelError',true,'SaveBStepInf',true);
solverFlow.NonLinearLoop();
printUtils.finalize()


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