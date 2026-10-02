function result = main(varargin)
% Dirichlet boundary condition enforcement (thesis section 3.3.2).
if nargin == 1 && isstruct(varargin{1}) && isfield(varargin{1}, 'SuiteRoot')
    result = runExperiment(varargin{1});
    return
end
suiteRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(suiteRoot, 'Utils'));
result = thesisRun(suiteRoot, 'Cases', "DirichletBoundaryConditionEnforcement", varargin{:});
end

function result = runExperiment(options)
% Migrated from Chapter_03/infsup/poissonApproxBC_clean.m.
result = struct();
%% poissonApproxBC_clean.m
% Poisson mesh-tying comparison on a 3D two-patch problem using GReS.
%
% This script runs all multiplier choices:
%   - standard
%   - dual
%   - P0
%
% and produces one plot for the primary variable u and one plot for the
% multiplier lambda along the slave-side interface diagonal.
%
% This version uses the standard GReS NonLinearImplicit solver.

%% User parameters
cfg.Nm      = 20;     % master elements in x/y directions
cfg.Nnm     = 20;     % slave/non-mortar elements in x/y directions
cfg.NzM     = 10;     % master elements in z direction
cfg.NzS     = 10;     % slave elements in z direction

cfg.multTypes = ["standard", "dual", "P0"];

cfg.variable      = "u";
cfg.masterSurface = 1;
cfg.slaveSurface  = 2;

cfg.stabilizationScale = 1e0;

cfg.tolDiagonal = 1e-5;
cfg.nAnalytical = 1000;

cfg.outputPrefix = "poisson_approxBC";

user = jsondecode(fileread('Input/config.json'));
for field = string(fieldnames(user))'
    cfg.(field) = user.(field);
end
if options.Smoke
    cfg.Nm = 4;
    cfg.Nnm = 4;
    cfg.NzM = 2;
    cfg.NzS = 2;
end

%% Manufactured solution
alpha = 1.0;

cfg.usol = @(x, y, z) sin(pi * x) .* sin(pi * y) .* (1 + alpha * (z - 0.5));

% Sign convention used by the current GReS Poisson setup in your scripts.
% If your local Poisson class is changed to -Delta u = f, flip the sign.
cfg.fsol = @(x, y, z) ...
  -2 * pi^2 .* sin(pi * x) .* sin(pi * y) .* (1 + alpha * (z - 0.5));

cfg.lambdasol = @(x, y) alpha * sin(pi * x) .* sin(pi * y);

%% Run all cases
setupPlotStyle();

res = cell(numel(cfg.multTypes), 1);

for k = 1:numel(cfg.multTypes)
    res{k} = runApproxCase(cfg, cfg.multTypes(k));
end

%% Plots
plotComparison(res, cfg, "u");
plotComparison(res, cfg, "lambda");

result.profiles = res;
result.configuration = cfg;

%% Local functions

end

function out = runApproxCase(cfg, multType)
fprintf("\nRunning approximate-BC case with multiplier type: %s\n", multType);

fixMasterTags = setdiff(1:6, cfg.masterSurface);
fixSlaveTags  = setdiff(1:6, cfg.slaveSurface);

gMaster = structuredMesh(cfg.Nm, cfg.Nm, cfg.NzM, [0 1], [0 1], [0.5 1]);
gSlave  = structuredMesh(cfg.Nnm, cfg.Nnm, cfg.NzS, [0 1], [0 1], [0 0.5]);

nListMaster = getBoundaryNodeList(gMaster, fixMasterTags);
nListSlave  = getBoundaryNodeList(gSlave, fixSlaveTags);

bcMaster = Boundaries(gMaster);
bcMaster.addBC('name', "fixMaster", ...
  'type', "dirichlet", ...
  'field', "node", ...
  'variable', cfg.variable, ...
  'entityListType', "bcList", ...
  'entityList', nListMaster);

UM = cfg.usol(gMaster.coordinates(nListMaster, 1), ...
  gMaster.coordinates(nListMaster, 2), ...
  gMaster.coordinates(nListMaster, 3));
bcMaster.addBCEvent("fixMaster", 'time', 0.0, 'value', UM);

bcSlave = Boundaries(gSlave);
bcSlave.addBC('name', "fixSlave", ...
  'type', "dirichlet", ...
  'field', "node", ...
  'variable', cfg.variable, ...
  'entityListType', "bcList", ...
  'entityList', nListSlave);

US = cfg.usol(gSlave.coordinates(nListSlave, 1), ...
  gSlave.coordinates(nListSlave, 2), ...
  gSlave.coordinates(nListSlave, 3));
bcSlave.addBCEvent("fixSlave", 'time', 0.0, 'value', US);

domainMaster = Discretizer('boundaries', bcMaster, 'grid', gMaster);
domainMaster.addPhysicsSolver('Poisson', 'gaussOrder', 4);

domainSlave = Discretizer('boundaries', bcSlave, 'grid', gSlave);
domainSlave.addPhysicsSolver('Poisson', 'gaussOrder', 4);

domains = [domainMaster; domainSlave];

domainMaster.getPhysicsSolver("Poisson").setAnalSolution(cfg.usol, cfg.fsol);
domainSlave.getPhysicsSolver("Poisson").setAnalSolution(cfg.usol, cfg.fsol);

domainMaster.initialize();
domainSlave.initialize();

interfInput = struct('masterDomain', 1, ...
  'slaveDomain', 2, ...
  'masterSurface', cfg.masterSurface, ...
  'slaveSurface', cfg.slaveSurface, ...
  'multiplierType', multType, ...
  'stabilizationScale', cfg.stabilizationScale);

interfaces = InterfaceSolver.add("MeshTying", domains, interfInput);
interface = interfaces{1};

outState = OutState('outputFile', fullfile('Output', "results_" + multType), 'printTimes', 1);

solver = NonLinearImplicit( ...
  'simulationparameters', SimulationParameters('Input/simParam.xml'), ...
  'output', outState, ...
  'domains', domains, ...
  'interface', interfaces);

solver.simulationLoop();

uSlave = domainSlave.getState.u;
lambda = interface.getState.multipliers;

[sU, uDiag] = extractPrimaryDiagonal(gSlave, domainSlave, uSlave, cfg.tolDiagonal);
[sL, lambdaDiag] = extractMultiplierDiagonal(interface, gSlave, bcSlave, lambda, ...
  multType, false, cfg.tolDiagonal);

out.multType = multType;
out.sU = sU;
out.u = uDiag;
out.sLambda = sL;
out.lambda = lambdaDiag;

fprintf("  diagonal u points      : %d\n", numel(out.sU));
fprintf("  diagonal lambda points : %d\n", numel(out.sLambda));
end

function [s, uDiag] = extractPrimaryDiagonal(g, domain, uFull, tol)
c = g.coordinates;
zInterface = max(c(:, 3));

ids = find(abs(c(:, 1) - c(:, 2)) <= tol & abs(c(:, 3) - zInterface) <= tol);

s = sqrt(c(ids, 1).^2 + c(ids, 2).^2);
[s, ord] = sort(s);
ids = ids(ord);

dofs = getLocalDoF(domain.dofm, 1, ids);
uDiag = uFull(dofs);
end

function [s, lambdaDiag] = extractMultiplierDiagonal(interface, gSlave, bcSlave, lambda, multType, modifyBoundaryMultipliers, tol)
gsint = interface.grids(MortarSide.slave);

if string(multType) == "P0"
    c = gsint.surfaces.center;
    ids = find(abs(c(:, 1) - c(:, 2)) <= tol);

    s = sqrt(c(ids, 1).^2 + c(ids, 2).^2);
    [s, ord] = sort(s);
    ids = ids(ord);

    lambdaDiag = lambda(ids);
    return
end

nodesGlobal = unique(gsint.surfaces.loc2glob(:), 'stable');
nodesGlobal = nodesGlobal(nodesGlobal > 0);

if modifyBoundaryMultipliers
    bcNodes = bcSlave.getTargetEntities("fixSlave");
    isBoundary = ismember(nodesGlobal, bcNodes);
    nodesGlobal = nodesGlobal(~isBoundary);
end

if numel(lambda) == numel(nodesGlobal)
    c = gSlave.coordinates(nodesGlobal, :);
else
    % Some MeshTying configurations store nodal multiplier locations
    % directly on the interface grid. Use this only as a fallback.
    c = gsint.coordinates;
    if size(c, 1) ~= numel(lambda)
        error("Cannot match nodal multipliers to coordinates: numel(lambda) = %d.", numel(lambda));
    end
end

ids = find(abs(c(:, 1) - c(:, 2)) <= tol);

s = sqrt(c(ids, 1).^2 + c(ids, 2).^2);
[s, ord] = sort(s);
ids = ids(ord);

lambdaDiag = lambda(ids);
end

function plotComparison(res, cfg, quantity)
[style, labels] = caseStyle();

xAna = linspace(0, 1, cfg.nAnalytical);
sAna = sqrt(xAna.^2 + xAna.^2);

fig = figure('Name', sprintf('%s comparison', quantity));
hold on;
grid on;
box on;

switch string(quantity)
    case "u"
        yAna = cfg.usol(xAna, xAna, 0.5);
        yLabel = '$u$';
        fileName = cfg.outputPrefix + "_u_diagonal.pdf";

    case "lambda"
        yAna = cfg.lambdasol(xAna, xAna);
        yLabel = '$\lambda$';
        fileName = cfg.outputPrefix + "_lambda_diagonal.pdf";

    otherwise
        error("Unknown quantity: %s", quantity);
end

plot(sAna, yAna, 'k-', 'LineWidth', 3, 'DisplayName', 'Analytical');

for k = 1:numel(res)
    rk = res{k};
    mt = char(rk.multType);

    switch string(quantity)
        case "u"
            x = rk.sU;
            y = rk.u;
        case "lambda"
            x = rk.sLambda;
            y = rk.lambda;
    end

    plot(x, y, ...
         'LineStyle', '-', ...
         'Marker', style.(mt).marker, ...
         'Color', style.(mt).color, ...
         'MarkerFaceColor', style.(mt).color, ...
         'MarkerSize', 6, ...
         'LineWidth', 0.8, ...
         'DisplayName', labels.(mt));
end

xlabel('coordinate (diagonal)', 'Interpreter', 'latex');
ylabel(yLabel, 'Interpreter', 'latex');
legend('Location', 'best', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
exportgraphics(fig, fullfile('Output', fileName), 'ContentType', 'vector');
end

function [style, labels] = caseStyle()
style.standard.color = [0.0000 0.4470 0.7410];
style.standard.marker = 'o';

style.dual.color = [0.4940 0.1840 0.5560];
style.dual.marker = 'd';

style.P0.color = [0.8500 0.3250 0.0980];
style.P0.marker = 's';

labels.standard = 'Standard';
labels.dual = 'Dual';
labels.P0 = 'P0';
end

function setupPlotStyle()
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 13);
end

function nodeList = getBoundaryNodeList(g, tags)
mask = ismember(g.surfaces.tag, tags);
nodeList = unique(g.surfaces.connectivity(mask, :));
nodeList = nodeList(nodeList > 0);
end
