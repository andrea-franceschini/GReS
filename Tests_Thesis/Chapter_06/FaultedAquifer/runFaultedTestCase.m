%% Faulted aquifer: flow first, contact mechanics second (thesis 6.3.3).
% Run initGReS before this script. All meshes are generated in memory by GReS.
caseRoot = fileparts(mfilename('fullpath'));
cd(caseRoot);
addpath(fullfile(caseRoot, 'Utils'));
assert(exist('structuredMesh', 'file') == 2, 'Run initGReS before this script.');

% Geometry and discretization: Figure 6.10, distances in metres.
config.length = 2000;
config.width = 1000;
config.thickness = 100;
config.faultPosition = 400;             % x at mid-depth
config.faultAngle = 40;                 % degrees from vertical
config.rockCells = [3 24 12];
config.soilCells = [12 32 22];
config.sandFraction = 0.60;
config.layerWarp = 10;                  % gentle horizon undulation, metres
config.faultMeshShift = 0.12;           % interior layer offset / thickness
config.wellDistance = 500;              % x distance from the fault
config.wellY = [0.4 0.6] * config.width;
config.wellDepthFraction = 0.30;        % centre of the sandy interval

% Flow data: years, kPa; reported pressure is positive drawdown.
config.wellTimes = [0 3 7 10];
config.wellDrawdown = [0 120 360 720];
config.outputTimes = 0:10;

% Mechanical data: MPa, metres; initial stress is compression-negative.
config.specificWeight = 0.021;         % MPa/m
config.K0 = 1 - sind(30);
config.frictionAngle = 30;
config.cohesion = 0.1;
config.stabilizationScale = 1;

% Replace the previous results; no numbered or dated run folders.
if isfolder('Output')
    rmdir('Output', 's');
end
mkdir('Output');
[mesh, meshInfo] = generateFaultedAquiferGrid(config);
plotFaultedAquiferGrid(mesh, meshInfo, config);
save('Output/mesh.mat', 'mesh', 'meshInfo', 'config', '-v7.3');

%% Stage 1: transient flow. The fault is hydraulically sealed.
[mesh, pressures, times] = runFlowSimulation(mesh, meshInfo, config);

%% Stage 2: prescribed drawdown drives contact mechanics on the same cells.
mechanical = runContactMechanicsSimulation(mesh, pressures, times, config);
save('Output/result.mat', 'mechanical', 'config', '-v7.3');
