function result = main(varargin)
% Stabilized piecewise constant multiplier (thesis section 3.2.3).
if nargin == 1 && isstruct(varargin{1}) && isfield(varargin{1}, 'SuiteRoot')
    result = runExperiment(varargin{1});
    return
end
suiteRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(suiteRoot, 'Utils'));
result = thesisRun(suiteRoot, 'Cases', "StabilizationScale", varargin{:});
end

function result = runExperiment(options)
% Migrated from Chapter_03/infsup/poissonApproxBC_P0_stabilizationSweep.m.
result = struct();
%% poissonApproxBC_P0_stabilizationSweep.m
% Poisson mesh-tying comparison for P0 multipliers with different
% stabilization scales, using the standard GReS NonLinearImplicit solver.
%
% The script runs one P0 case for each user-selected stabilization scale and
% produces two plots along the slave-side interface diagonal:
%   1) primary variable u;
%   2) multiplier lambda.
%
% Analytical solution: continuous black line.
% Numerical solutions: filled markers only, no connecting line.

%% User parameters
cfg.Nm      = 20;     % master elements in x/y directions
cfg.Nnm     = 20;     % slave/non-mortar elements in x/y directions
cfg.NzM     = 10;     % master elements in z direction
cfg.NzS     = 10;     % slave elements in z direction

cfg.multType = "P0";

% User-selected stabilization scales.
cfg.stabilizationScales = [0.0, 1e0, 1e1, 1e2];

cfg.variable      = "u";
cfg.masterSurface = 1;
cfg.slaveSurface  = 2;

cfg.tolDiagonal = 1e-5;
cfg.nAnalytical = 1000;

cfg.outputPrefix = "poisson_approxBC_P0_stabilizationSweep";

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

%% Run all stabilization cases
setupPlotStyle();

res = cell(numel(cfg.stabilizationScales), 1);
for k = 1:numel(cfg.stabilizationScales)
    try
        res{k} = runApproxP0Case(cfg, cfg.stabilizationScales(k));
        res{k}.success = all(isfinite([res{k}.u(:); res{k}.lambda(:)]));
        res{k}.message = "";
    catch ME
        res{k} = struct('stabilizationScale', cfg.stabilizationScales(k), ...
         'sU', [], 'u', [], 'sLambda', [], 'lambda', [], 'success', false, ...
         'message', string(getReport(ME, 'extended', 'hyperlinks', 'off')));
        warning('Thesis:StabilizationFailure', 'Scale %g failed: %s', cfg.stabilizationScales(k), ME.message);
    end
end

%% Plots
plotScaleComparison(res, cfg, "u");
plotScaleComparison(res, cfg, "lambda");

result.failedConfigurations = sum(cellfun(@(r) ~r.success, res));
result.profiles = res;
result.configuration = cfg;

%% Local functions

end

function out = runApproxP0Case(cfg, stabilizationScale)
fprintf("\nRunning approximate-BC P0 case with stabilization scale: %.6e\n", stabilizationScale);

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
                     'multiplierType', cfg.multType, ...
                     'stabilizationScale', stabilizationScale);

interfaces = InterfaceSolver.add("MeshTying", domains, interfInput);
interface = interfaces{1};

outState = OutState('outputFile', fullfile('Output', "results_P0_stab_" + scaleTag(stabilizationScale)), ...
                    'printTimes', 1);

solver = NonLinearImplicit( ...
    'simulationparameters', SimulationParameters('Input/simParam.xml'), ...
    'output', outState, ...
    'domains', domains, ...
    'interface', interfaces);

solver.simulationLoop();

uSlave = domainSlave.getState.u;
lambda = interface.getState.multipliers;

[sU, uDiag] = extractPrimaryDiagonal(gSlave, domainSlave, uSlave, cfg.tolDiagonal);
[sL, lambdaDiag] = extractP0MultiplierDiagonal(interface, lambda, cfg.tolDiagonal);

out.multType = cfg.multType;
out.stabilizationScale = stabilizationScale;
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

function [s, lambdaDiag] = extractP0MultiplierDiagonal(interface, lambda, tol)
gsint = interface.grids(MortarSide.slave);

c = gsint.surfaces.center;
ids = find(abs(c(:, 1) - c(:, 2)) <= tol);

s = sqrt(c(ids, 1).^2 + c(ids, 2).^2);
[s, ord] = sort(s);
ids = ids(ord);

lambdaDiag = lambda(ids);
end

function plotScaleComparison(res, cfg, quantity)
xAna = linspace(0, 1, cfg.nAnalytical);
sAna = sqrt(xAna.^2 + xAna.^2);

colors = lines(numel(res));
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', 'x', '+'};

fig = figure('Name', sprintf('P0 stabilization scale comparison: %s', quantity));
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

plot(sAna, yAna, 'k-', 'LineWidth', 1.4, 'DisplayName', 'Analytical');

for k = 1:numel(res)
    rk = res{k};

    switch string(quantity)
        case "u"
            x = rk.sU;
            y = rk.u;
        case "lambda"
            x = rk.sLambda;
            y = rk.lambda;
    end

    marker = markers{mod(k - 1, numel(markers)) + 1};
    label = sprintf('$\\alpha = %g$', rk.stabilizationScale);

    plot(x, y, ...
         'LineStyle', '-', ...
         'Marker', marker, ...
         'Color', colors(k, :), ...
         'MarkerFaceColor', colors(k, :), ...
         'MarkerSize', 5.5, ...
         'LineWidth', 0.8, ...
         'DisplayName', label);
end

xlabel('coordinate (diagonal)', 'Interpreter', 'latex');
ylabel(yLabel, 'Interpreter', 'latex');
legend('Location', 'best', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
ylim([-0.1, 1.2]);
exportgraphics(fig, fullfile('Output', fileName), 'ContentType', 'vector');
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

function tag = scaleTag(scale)
tag = string(sprintf('%.3e', scale));
tag = replace(tag, "+", "");
tag = replace(tag, "-", "m");
tag = replace(tag, ".", "p");
end
