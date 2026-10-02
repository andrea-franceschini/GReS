function result = main(varargin)
% Fluid withdrawal from a deep aquifer (thesis section 6.3.1).
if nargin == 1 && isstruct(varargin{1}) && isfield(varargin{1}, 'SuiteRoot')
    result = runExperiment(varargin{1});
    return
end
suiteRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(suiteRoot, 'Utils'));
result = thesisRun(suiteRoot, 'Cases', "FluidWithdrawalFromADeepAquifer", varargin{:});
end

function result = runExperiment(options)
% Migrated from Chapter_06/FancyMultiDomain/deepAquifer_3.m.
result = struct();

% mesh processing
fprintf('Processing meshes...\n');

meanDepth     = 700;   % mean depth of upper surface [m], positive downward
meanThickness = 150;    % mean aquifer thickness [m]

surfaceVar    = 0.35;   % variability of upper surface: 0 = flat, 1 = full
thicknessVar  = 0.30;   % variability of thickness: 0 = constant, 1 = full

dipX = -150;            % total elevation change along x [m]
dipY =  110;            % total elevation change along y [m]

minThickness = 60;     % minimum allowed thickness [m]

nx = 20;
if options.Smoke
    nx = 8;
end
ny = nx;

x = linspace(-1e3, 1e3, nx);
y = linspace(-1e3, 1e3, ny);

[X, Y] = ndgrid(x, y);

xN = X / max(x);
yN = Y / max(y);

xp = 6 * (xN - 0.5);
yp = 6 * (yN - 0.5);

%% Upper-surface shape

upperVariation = ...
      dipX * (xN - 0.5) ...
    + dipY * (yN - 0.5) ...
    + 120 * sin(2 * pi * xN) .* cos(pi * yN) ...
    + 20 * peaks(xp, yp);

% Remove its mean so mean(zUpper) is exactly -meanDepth
upperVariation = upperVariation - mean(upperVariation, 'all');

zUpper = -meanDepth + surfaceVar * upperVariation;

%% Thickness field

thicknessVariation = ...
      70 * sin(pi * xN).^2 ...
    + 40 * cos(2 * pi * yN) ...
    + 25 * exp(-((xp - 1).^2 + (yp + 0.5).^2));

% Remove the mean so the requested mean thickness is preserved
thicknessVariation = ...
    thicknessVariation - mean(thicknessVariation, 'all');

thickness = meanThickness ...
          + thicknessVar * thicknessVariation;

thickness = max(thickness, minThickness);

%% Lower surface

zLower = zUpper - thickness;
% Use the native structured generator and map its vertical layers between
% the same analytical horizons as the original MRST construction.
aquiferMesh = structuredMesh(x, y, linspace(0, 1, 8));
c = aquiferMesh.coordinates;
zlo = interp2(X', Y', zLower', c(:, 1), c(:, 2), 'linear');
zhi = interp2(X', Y', zUpper', c(:, 1), c(:, 2), 'linear');
aquiferMesh.coordinates(:, 3) = zlo + c(:, 3) .* (zhi - zlo);

%%

burdenMesh = structuredMesh([4 15 4], [4 15 4], [3 15 4], [-2e3, -1.01e3, 1.01e3, 2e3], [-2e3, -1.01e3, 1.01e3, 2e3], [-2e3 -1e3 -6e2 0]);

burdenMesh = burdenMesh.gridDifference(aquiferMesh, 0.0);

aquiferMesh.processGeometry;
burdenMesh.processGeometry;

fprintf('Done Processing meshes...\n');

fprintf('Total nodes: %i \n', aquiferMesh.nNodes + burdenMesh.nNodes);
%%
% materials

matBurden = Materials();
matBurden.addSolid('name', "rock", 'cellTags', 1);
matBurden.addConstitutiveLaw("rock", "Elastic", 'youngModulus', 5e4, 'poissonRatio', 0.25);

matAquifer = Materials();
matAquifer.addFluid('dynamicViscosity', 3.1608e-14, 'specificWeight', 0.0, 'compressibility', 4.4e-7);
matAquifer.addSolid('name', "sand", 'cellTags', 1);
matAquifer.addPorousRock("sand", "permeability", 1e-13, "porosity", 0.25);
matAquifer.addConstitutiveLaw("sand", "Elastic", 'youngModulus', 6e4, 'poissonRatio', 0.25);

% boundary conditions

aquiferCenter = mean(aquiferMesh.coordinates, 1);
tol = 60;
wells = abs(aquiferMesh.cells.center(:, 1) - aquiferCenter(1)) < tol & ...
        abs(aquiferMesh.cells.center(:, 2) - aquiferCenter(2)) < tol;

fprintf('Pressure prescribed to %i aquifer cells\n', sum(wells));

c = burdenMesh.coordinates;
cMin = min(burdenMesh.coordinates, [], 1);
cMax = max(burdenMesh.coordinates, [], 1);

mBot = find(abs(c(:, 3) - cMin(3)) < 1e-3);
mLatX = find(abs(c(:, 1) - cMax(1)) < 1e-3 | abs(c(:, 1) - cMin(1)) < 1e-3);
mLatY = find(abs(c(:, 2) - cMax(2)) < 1e-3 | abs(c(:, 2) - cMin(2)) < 1e-3);

c = aquiferMesh.cells.center;
cMin = min(aquiferMesh.cells.center, [], 1);
cMax = max(aquiferMesh.cells.center, [], 1);
boundAquiferCells = find(abs(c(:, 2) - cMin(2)) < 100 | abs(c(:, 2) - cMax(2)) < 1e2 | abs(c(:, 1) - cMin(1)) < 1e2 | abs(c(:, 1) - cMax(1)) < 1e2);

bcAquifer = Boundaries(aquiferMesh);

bcAquifer.addBC('name', "pressureWells", ...
        'type', "source", ...
        'field', "cell", ...
        'entityListType', "bcList", ...
        'entityList', find(wells), ...
        'variable', "pressure");
bcAquifer.addBCEvent("pressureWells", 'time', 0.0, 'value', -0.1);

bcAquifer.addBC('name', "pressureBound", ...
        'type', "dirichlet", ...
        'field', "cell", ...
        'entityListType', "bcList", ...
        'entityList', boundAquiferCells, ...
        'variable', "pressure");
bcAquifer.addBCEvent("pressureBound", 'time', 0.0, 'value', 0.0);

bcBurden = Boundaries(burdenMesh);

bcBurden.addBC('name', "fixX", ...
        'type', "dirichlet", ...
        'field', "node", ...
        'entityListType', "bcList", ...
        'entityList', mLatX, ...
        'components', "x", ...
        'variable', "displacements");

bcBurden.addBCEvent("fixX", 'time', 0.0, 'value', 0.0);

bcBurden.addBC('name', "fixY", ...
        'type', "dirichlet", ...
        'field', "node", ...
        'entityListType', "bcList", ...
        'entityList', mLatY, ...
        'components', "y", ...
        'variable', "displacements");

bcBurden.addBCEvent("fixY", 'time', 0.0, 'value', 0.0);

bcBurden.addBC('name', "fixBot", ...
        'type', "dirichlet", ...
        'field', "node", ...
        'entityListType', "bcList", ...
        'entityList', mBot, ...
        'components', "z", ...
        'variable', "displacements");

bcBurden.addBCEvent("fixBot", 'time', 0.0, 'value', 0.0);

% domains and interface

fprintf('Processing domains and interface...\n');

domBurden = Discretizer('grid', burdenMesh, 'boundaries', bcBurden, 'materials', matBurden);
domAquifer = Discretizer('grid', aquiferMesh, 'boundaries', bcAquifer, 'materials', matAquifer);

domBurden.addPhysicsSolver('Poromechanics');
domAquifer.addPhysicsSolver('BiotFullyCoupled');
domains = [domBurden, domAquifer];

interf.masterDomain = 1;
interf.slaveDomain = 2;
interf.scale = 0.25;
interf.Quadrature.type = "RBFquadrature";
interface = InterfaceSolver.add('MeshTying', domains, interf);
fprintf('Done Processing domains and interface...\n');
g = interface{1}.grids;
fprintf('Found %i master elements and %i slave elements \n', g(2).surfaces.num, g(1).surfaces.num);

printUtils = OutState('outputFile', "Output/deepAquifer", 'printTimes', 0:0.1:1);

simparams = SimulationParameters("Input/simParam.xml");

solver = NonLinearImplicit('simulationparameters', simparams, 'domains', domains, 'output', printUtils, 'interface', interface);
solver.simulationLoop();

result.wellCells = find(wells);
result.finalTime = 1;
%%

end
