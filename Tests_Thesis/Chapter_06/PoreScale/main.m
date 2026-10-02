%% Flow problem

fprintf('Flow problem... \n\n');
fluidMesh = Grid();
fluidMesh.importMesh('Mesh/background.vtk');

% get surfaces
% 1 - internal grain surfaces
% 2 - bottom surfaces
% 3 - top surfaces
% fluidMesh = setSurfacesFluid(fluidMesh);

input = struct('Start', 0.0, 'End', 1.0, 'DtInit', 1.0, 'DtMax', 1.0, 'DtMin', 1.0, 'AbsoluteTolerance', 1e-6, 'RelativeTolerance', 1e-8);
simparams = SimulationParameters(input);

fluidMesh.processGeometry();

% spot cells that do not have any neighboring cells (no internal faces)
f = fluidMesh.faces;
intFaces = all(f.neighbors, 2);
validCells = unique(f.neighbors(intFaces, :));
id = false(fluidMesh.cells.num, 1);
id(validCells) = true;
fluidMesh = getCellGrid(fluidMesh, id);
fluidMesh.processGeometry();

mat = Materials();
mat.addFluid('dynamicViscosity', 1e-6, 'specificWeight', 0.0, 'compressibility', 0.0);
mat.addSolid('name', "sand", 'cellTags', 1);
mat.addPorousRock("sand", "permeability", 1e-10, "porosity", 0.375);

% get bottom and top surfaces
botBox = [0.1499, 0.8501, 0.1499, 0.8501, 0.1499, 0.1501]; % box inscribing the bottom surface
topBox = [0.1499, 0.8501, 0.1499, 0.8501, 0.8499, 0.8501]; % box inscribing the bottom surface

bc = Boundaries(fluidMesh);
bc.addBC('name', "bottomPressure", ...
        'type', "dirichlet", ...
        'field', "surface", ...
        'entityListType', "box", ...
        'entityList', botBox, ...
        'variable', "pressure");

% top surfaces have tag 2 in structuredMesh
bc.addBC('name', "topPressure", ...
        'type', "dirichlet", ...
        'field', "surface", ...
        'entityListType', "box", ...
        'entityList', topBox, ...
        'variable', "pressure");

% bottom pressure - single event -> constant
bc.addBCEvent("bottomPressure", 'time', 0.0, 'value', 0.0);
bc.addBCEvent("topPressure", 'time', 0.0, 'value', 1e3);

% discretizer
domFluid = Discretizer('Boundaries', bc, ...
                       'Materials',  mat, ...
                       'Grid',       fluidMesh);

domFluid.addPhysicsSolver("SinglePhaseFlowFVTPFA");

out = OutState('printTimes', 0:1:100, 'outputFile', "Output/fluidPore", 'matFileName', "Output/fluidPore", 'vtkFormat', "ascii");

solver = NonLinearImplicit('simulationparameters', simparams, 'domains', domFluid, 'output', out);
gresLog().setVerbosity(2);
solver.simulationLoop();

%% interpolate pressure from fluid mesh to grain surfaces

fprintf('Interpolating pressure from void space to grain surfaces... \n\n');
% import grain mesh

grainMesh = Grid();
grainMesh.importMesh('Mesh/many_grains.vtk');
grainMesh.cells.tag = grainMesh.cells.tag + 1;
grainMesh.cells.nTag = 1;
grainMesh.processGeometry;

if isfile("interpMat.mat")

    interpMat = load("interpMat.mat");
    E = interpMat.E;
    slaveNodes = interpMat.slaveNodes;
    masterNoddes = interpMat.masterNodes;

else

    % compute the cross-grid interpolation operator

    % define the mortar interface
    domFluidNode = Discretizer('Grid', fluidMesh);
    domFluidNode.addPhysicsSolver('SinglePhaseFlowFEM');

    % create domains
    domGrain = Discretizer('Grid', grainMesh);
    domGrain.addPhysicsSolver('SinglePhaseFlowFEM');

    interfIn = struct('masterDomain', 1, 'slaveDomain', 2, 'multiplierType', "dual");

    tic;
    fprintf('Processing pore scale interface. Could take a lot... \n');
    interf = InterfaceSolver.add("MeshTying", [domFluidNode, domGrain], interfIn);
    t = toc;
    fprintf('Done Processing pore scale interface in %d s \n', t);

    interf{1}.computeConstraintMatrices();

    slaveNodes = interf{1, 1}.grids(MortarSide.slave).surfaces.loc2glob;
    masterNodes = interf{1, 1}.grids(MortarSide.master).surfaces.loc2glob;
    D = interf{1}.D(:, slaveNodes);
    M = interf{1}.M(:, masterNodes);

    E = D \ M;

    save('interpMat.mat', "E", '-mat');
    save('interpMat.mat', "masterNodes", '-mat', '-append');
    save('interpMat.mat', "slaveNodes", '-mat', '-append');

end

pressCell = domFluid.getState().pressure;

pressNode = zeros(fluidMesh.nNodes, 1);
volNode = zeros(fluidMesh.nNodes, 1);

h = Hexahedron(fluidMesh);

topol = fluidMesh.getCellNodes();

for i = 1:fluidMesh.cells.num
    nId = topol(i, :);
    volNod = h.getNodeInfluence(i);
    pressNode(nId) = pressNode(nId) + pressCell(i) * volNod;
    volNode(nId) = volNode(nId) + volNod;
end

areaNode = zeros(grainMesh.nNodes, 1);
t = Triangle(grainMesh);
topol = grainMesh.getSurfNodes();

for i = 1:grainMesh.surfaces.num
    nId = topol(i, :);
    areaNod = t.getNodeInfluence(i);
    areaNode(nId) = areaNode(nId) + areaNod;
end

pressNode = pressNode ./ volNode;

% interpolate pressure from flow mesh to grain surfaces
pressGrain = E * pressNode(masterNodes);

pressGrain(pressGrain > 10) = 10;
pressGrain(pressGrain < 0) = 0;

% avarage node pressure to get surface pressure
grainSurfTopol = grainMesh.getSurfNodes();

pGrainPlot = zeros(grainMesh.nNodes, 1);
pGrainPlot(slaveNodes) = pressGrain;

% get pressure on the surfaces of the grains
pressSurf = pGrainPlot(grainSurfTopol);
pressSurf = sum(pressSurf, 2) / 3;

%% setup mechanical problem

fprintf('Mechanical simulation \n\n');

% Find pairs of grain connectivities
cs = ContactSearching(grainMesh, grainMesh, 'scale', 1e-3);
elemConn = cs.getElementConnectivity();

[m, s] = find(elemConn);

% Identify connected grain pairs
cells = grainMesh.cells;
faces = grainMesh.faces;

% Cell adjacent to each detected surface
cM = sum(faces.neighbors(grainMesh.surfaces.faceId(m)), 2);
cS = sum(faces.neighbors(grainMesh.surfaces.faceId(s)), 2);

% Original grain tags
grainM = cells.tag(cM);
grainS = cells.tag(cS);

% Remove self-contact

isSelf = grainM == grainS;

grainM(isSelf) = [];
grainS(isSelf) = [];

% Remove duplicated grain pairs

% Represent every pair as (g1,g2), with g1 < g2
g1 = min(grainM, grainS);
g2 = max(grainM, grainS);

grainPairs = unique([g1, g2], 'rows');

% Identify mortar-active grains and renumber surface tags

activeGrains = unique(grainPairs(:));
nActiveGrains = numel(activeGrains);

% Last tag is reserved for surfaces belonging to grains that do not
% participate in any mortar interface
unconnectedTag = nActiveGrains + 1;

% Map:
%   original active grain tag -> 1,...,nActiveGrains
%   unconnected grain         -> nActiveGrains + 1
grainToSurfaceTag = unconnectedTag * ones(cells.nTag, 1);
grainToSurfaceTag(activeGrains) = 1:nActiveGrains;

% Find grain owning each surface
surfaceCells = sum( ...
    faces.neighbors(grainMesh.surfaces.faceId), ...
    2);

surfaceGrains = cells.tag(surfaceCells);

% Assign new surface tags
grainMesh.surfaces.tag = grainToSurfaceTag(surfaceGrains);
grainMesh.surfaces.nTag = unconnectedTag;

% Build mortar pairs using the NEW surface tags

% Convert grain-pair numbering to the new surface-tag numbering
surfacePairs = [ ...
    grainToSurfaceTag(grainPairs(:, 1)), ...
    grainToSurfaceTag(grainPairs(:, 2))];

% pairs{k} contains all master surfaces associated with slave surface k
pairsAll = cell(nActiveGrains, 1);

for k = 1:size(surfacePairs, 1)

    slaveTag  = surfacePairs(k, 1);
    masterTag = surfacePairs(k, 2);

    pairsAll{slaveTag}(end + 1) = masterTag;

end

% Keep only slave surfaces that actually have a mortar pair
slaveSurface = find(~cellfun(@isempty, pairsAll));

% Cell array now contains no empty entries
pairs = pairsAll(slaveSurface);

% Mechanical model

matGrain = Materials();

matGrain.addSolid( ...
    'name', "sand", ...
    'cellTags', 1:cells.nTag);

matGrain.addConstitutiveLaw( ...
    "sand", "Elastic", ...
    'youngModulus', 2e6, ...
    'poissonRatio', 0.25);

% Mechanical boundary conditions

bcGrain = Boundaries(grainMesh);

% Find external surfaces
cmin = 0.15;
cmax = 0.85;
tol  = 1e-3;

bot = abs(grainMesh.surfaces.center(:, 3) - cmin) < tol;
top = abs(grainMesh.surfaces.center(:, 3) - cmax) < tol;

south = abs(grainMesh.surfaces.center(:, 2) - cmin) < tol;
north = abs(grainMesh.surfaces.center(:, 2) - cmax) < tol;

west = abs(grainMesh.surfaces.center(:, 1) - cmin) < tol;
east = abs(grainMesh.surfaces.center(:, 1) - cmax) < tol;

surfs = [bot, top, south, north, west, east];

names = ["bot", "top", "south", "north", "west", "east"];

% Normal displacement component for each boundary
dir = [3, 3, 2, 2, 1, 1];

sBound = find(any(surfs, 2));

% Remaining surfaces receive pressure traction
sPress = setdiff((1:grainMesh.surfaces.num)', sBound);

% Fix normal displacement on external boundaries

for i = 1:6

    bcGrain.addBC( ...
        'name', names(i), ...
        'type', "dirichlet", ...
        'field', "surface", ...
        'entityListType', "bcList", ...
        'entityList', find(surfs(:, i)), ...
        'variable', "displacements", ...
        'components', dir(i));

    bcGrain.addBCEvent( ...
        names(i), ...
        'time', 0.0, ...
        'value', 0.0);

end

% Surface traction: t = p n

t = pressSurf(sPress) .* grainMesh.surfaces.normal(sPress, :);

tractionNames = ["tX", "tY", "tZ"];

for i = 1:3

    bcGrain.addBC( ...
        'name', tractionNames(i), ...
        'type', "neumann", ...
        'field', "surface", ...
        'entityListType', "bcList", ...
        'entityList', sPress, ...
        'variable', "displacements", ...
        'components', i);

    bcGrain.addBCEvent( ...
        tractionNames(i), ...
        'time', 0.0, ...
        'value', t(:, i));

end

% Mechanical domain

domGrain = Discretizer( ...
    'Boundaries', bcGrain, ...
    'Materials',  matGrain, ...
    'Grid',       grainMesh);

domGrain.addPhysicsSolver("Poromechanics");

% Add mortar interfaces

if ~isempty(pairs)

    interfaces = InterfaceSolver.add( ...
        "MeshTying", ...
        [domGrain, domGrain], ...
        struct('masterSurface', pairs{1}, 'slaveSurface', slaveSurface(1)));

    for k = 2:numel(pairs)

        interfaces = InterfaceSolver.add( ...
            "MeshTying", ...
            [domGrain, domGrain], ...
            interfaces, ...
             struct('masterSurface', pairs{k}, 'slaveSurface', slaveSurface(k)));

    end

end

% Output and solver
grainMesh.processGeometry;

out = OutState( ...
    'printTimes', 1, ...
    'outputFile', "Output/mechPore", ...
    'matFileName', "Output/mechPore", ...
    'vtkFormat', "ascii");

input = struct('Start', 0.0, 'End', 1.0, 'DtInit', 1.0, 'DtMax', 1.0, 'DtMin', 1.0, 'AbsoluteTolerance', 1e-4, 'RelativeTolerance', 1e-3);
simparams = SimulationParameters(input);

solver = NonLinearImplicit( ...
    'simulationparameters', simparams, ...
    'domains', domGrain, ...
    'output', out, ...
    'interface', interfaces);

gresLog().setVerbosity(2);

solver.simulationLoop();
