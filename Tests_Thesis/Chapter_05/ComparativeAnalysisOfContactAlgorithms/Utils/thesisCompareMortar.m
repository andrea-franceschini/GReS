function result = thesisCompareMortar(options)
% Migrated from Chapter_05/StickSlipThesis/Mortar_clean/Mortar/runComparison.m.
result = struct();

gresLog().setVerbosity(1);
%% Mortar stick-slip-open comparison
% The script creates two comparisons for both loading directions:
%   1. smooth, smooth-regularized and semi-smooth formulations at one
%      fixed augmentation scaling;
%   2. smooth-regularized and semi-smooth formulations while that scaling
%      changes.
% Before-bounce and after-bounce results are shown in separate panels.

%% User settings
% These values are dimensionless scale factors alpha. The dimensional
% mortar augmentation used by the formulation is c = alpha*E/h.
config.fixedAugmentation = 1e0;
config.augmentationSweep = [1e-2, 1e-1, 1e0, 1e1, 1e2];
% config.augmentationSweep = [1e0];
% h used by the existing mortar benchmark (the previous E/2 scaling).
config.characteristicLength = 1;
config.sweepMethodIds = ["smoothregularized", "semismooth"];
config.directions = ["vertical", "horizontal"];
config.scenarios = ["pre_bounce", "post_bounce"];
config.scenarioLabels = ["before bounce", "after bounce"];
config.endTimes = [11, 16];

% The first installed class in each list is used. The candidate lists keep
% the script usable with the naming variants found in GReS development forks.
config.methods = struct( ...
                        'id', {"smooth", "smoothregularized", "semismooth"}, ...
                        'solverCandidates', {{"SolidMechanicsContact"}, {"SolidMechanicsContactNew"}, ...
                       {"SolidMechanicsContactAugmented"}}, ...
                        'label', {"smooth", "smooth regularized", "semi-smooth"}, ...
                        'color', {[0.000, 0.447, 0.741], [0.929, 0.694, 0.125], [0.494, 0.184, 0.556]});

config.inputFile = fullfile("Input", "StickSlipOpen.xml");
config.materialFile = fullfile("Input", "materials.xml");
config.outputDir = ".";
config.historyDir = fullfile(config.outputDir, "Histories");
config.figureDir = fullfile(config.outputDir, "Figures");
config.resultsFile = fullfile(config.outputDir, "comparison_results.mat");

ensureDirectory(config.outputDir);
ensureDirectory(config.historyDir);
ensureDirectory(config.figureDir);

for iMethod = 1:numel(config.methods)
    config.methods(iMethod).solverName = resolveFormulation(config.methods(iMethod));
end

records = initializeRecords();
if isfile(config.resultsFile)
    saved = load(config.resultsFile, 'records');
    if isfield(saved, 'records')
        records = saved.records;
    end
end

%% Fixed-c formulation comparison
for iDirection = 1:numel(config.directions)
    direction = config.directions(iDirection);
    for iScenario = 1:numel(config.scenarios)
        scenario = config.scenarios(iScenario);
        endTime = config.endTimes(iScenario);
        for iMethod = 1:numel(config.methods)
            method = config.methods(iMethod);
            [records, wasRun] = runAndRecord(records, config, direction, ...
                                             scenario, endTime, method, config.fixedAugmentation);
            if wasRun
                save(config.resultsFile, 'records', 'config', '-v7.3');
            end
        end
    end
end

%% Smooth-regularized and semi-smooth augmentation sweeps
methodIds = string({config.methods.id});
sweepMethods = config.methods(ismember(methodIds, ...
                                       config.sweepMethodIds));

for iDirection = 1:numel(config.directions)
    direction = config.directions(iDirection);

    for iScenario = 1:numel(config.scenarios)
        scenario = config.scenarios(iScenario);
        endTime = config.endTimes(iScenario);

        for iMethod = 1:numel(sweepMethods)
            method = sweepMethods(iMethod);

            for augmentation = config.augmentationSweep
                [records, wasRun] = runAndRecord(records, config, direction, ...
                                                 scenario, endTime, method, augmentation);

                if wasRun
                    save(config.resultsFile, 'records', 'config', '-v7.3');
                end
            end
        end
    end
end

%% Figures
for direction = config.directions
    plotFixedComparison(records, config, direction);
    plotAugmentationSweep(records, config, direction, sweepMethods);
end

fprintf('\nMortar comparison complete.\n');
fprintf('Results: %s\n', config.resultsFile);
fprintf('Figures: %s\n', config.figureDir);

%% Local functions

result.records = records;
result.failedConfigurations = sum(~[records.success]);
end

function records = initializeRecords()
records = struct( ...
                 'direction', {}, 'scenario', {}, 'method', {}, 'solverName', {}, ...
                 'augmentation', {}, 'iterations', {}, 'success', {}, 'message', {}, ...
                 'historyFile', {});
end

function solverName = resolveFormulation(method)
solverName = "";
candidates = string(method.solverCandidates);
for candidate = candidates
    if exist(candidate, 'class') == 8
        solverName = candidate;
        return
    end
end
warning('Mortar:MissingFormulation', ...
        ['No installed class was found for "%s". The run will be recorded as ', ...
             'failed. Edit config.methods.solverCandidates for your GReS fork.'], ...
        method.label);
solverName = candidates(1);
end

function [records, wasRun] = runAndRecord(records, config, direction, ...
                                          scenario, endTime, method, augmentation)

idx = findRecord(records, direction, scenario, method.id, augmentation);
% if ~isempty(idx) && records(idx).success
%  fprintf('Using cached result: %s, %s, %s, c/(E/h)=%g\n', ...
%    direction, scenario, method.label, augmentation);
%  wasRun = false;
%  return
% end
wasRun = true;

fprintf('\n%s | %s | %s | c/(E/h) = %.4e\n', ...
        upper(direction), scenario, method.label, augmentation);

tag = sprintf('%s_%s_%s_c_%s', direction, scenario, method.id, ...
              scientificTag(augmentation));
historyBase = fullfile(config.historyDir, "Mortar_" + tag);

entry.direction = direction;
entry.scenario = scenario;
entry.method = method.id;
entry.solverName = method.solverName;
entry.augmentation = augmentation;
entry.iterations = NaN;
entry.success = false;
entry.message = "";
entry.historyFile = historyBase + ".mat";

try
    [solver, centers] = runCase(config, direction, endTime, ...
                                method.solverName, augmentation, historyBase);
    entry.iterations = solver.totIter;
    entry.success = true;
    appendPlotStepMetadata(entry.historyFile, centers, direction, ...
                           scenario, method.label, augmentation, "Mortar");
    fprintf('iterations = %g\n', entry.iterations);
catch ME
    entry.success = false;
    entry.message = string(ME.getReport('extended', 'hyperlinks', 'off'));
    fprintf(2, 'FAILED: %s\n', ME.message);
end

if isempty(idx)
    records(end + 1) = entry; %#ok<AGROW>
else
    records(idx) = entry;
end
wasRun = true;
end

function idx = findRecord(records, direction, scenario, method, augmentation)
idx = [];
if isempty(records)
    return
end
same = [records.direction] == direction & ...
  [records.scenario] == scenario & ...
  [records.method] == method & ...
  abs([records.augmentation] - augmentation) <= ...
  10 * eps(max(1, abs(augmentation)));
idx = find(same, 1, 'last');
end

function [solver, centers] = runCase(config, direction, endTime, ...
                                     formulation, augmentation, historyBase)

if exist(formulation, 'class') ~= 8
    error('Mortar:MissingFormulation', ...
          ['The GReS formulation class "%s" is not on the MATLAB path. ', ...
                 'Update config.methods.solverCandidates if your fork uses another name.'], ...
          formulation);
end

params = readInput(config.inputFile);
params.SimulationParameters.End = endTime;
simParam = SimulationParameters(params.SimulationParameters);

X = 5;
Y = 10;
Z = 15;
nxLeft = 3;
nyLeft = 8;
nzLeft = 12;
nxRight = 3;
nyRight = 8;
nzRight = 12;
gridLeft = structuredMesh(nxLeft, nyLeft, nzLeft, ...
                          [0, 0.5 * X], [0, Y], [0, Z]);
gridRight = structuredMesh(nxRight, nyRight, nzRight, ...
                           [0.5 * X, X], [0, Y], [0, Z]);

assert(mod(nyLeft, 2) == 0 && mod(nyRight, 2) == 0, ...
       'The number of elements along the y axis must be even.');

materialsLeft = Materials(config.materialFile);
materialsRight = Materials(config.materialFile);

E = materialsLeft.getConstitutiveLaw(1).E;
augmentation = augmentation * E / config.characteristicLength;

switch direction
    case "vertical"
        [bcLeft, bcRight] = setVerticalBC(Y, gridLeft, gridRight);
    case "horizontal"
        [bcLeft, bcRight] = setHorizontalBC(Y, gridLeft, gridRight);
    otherwise
        error('Mortar:LoadingDirection', 'Unknown loading direction: %s', direction);
end

domainLeft = Discretizer( ...
                         'boundaries', bcLeft, 'materials', materialsLeft, 'grid', gridLeft);
domainRight = Discretizer( ...
                          'boundaries', bcRight, 'materials', materialsRight, 'grid', gridRight);
domainLeft.addPhysicsSolver('Poromechanics');
domainRight.addPhysicsSolver('Poromechanics');
domains = [domainLeft; domainRight];

interfaceNames = fieldnames(params.Interface);
if numel(interfaceNames) ~= 1
    error('Mortar:InputInterface', ...
          'The comparison input must define exactly one interface template.');
end
interfaceInput = params.Interface.(interfaceNames{1});
interfaceInput.augmentationParameter = augmentation;
interfaceInput.augmentation = augmentation;
interfaceInput.augmentationNormal = augmentation;
interfaceInput.augmentationTangential = augmentation;
interfaceInput.penaltyParameter = augmentation;
interfaceInput.penaltyNormal = augmentation;
interfaceInput.penaltyTangential = augmentation;
params.Interface = struct();
params.Interface.(char(formulation)) = interfaceInput;
interfaces = InterfaceSolver.addInterfaces(domains, params.Interface);

stressLeft = getState(domainLeft, "stress");
stressLeft(:, 1) = -1.0;
setState(domainLeft, stressLeft, "stress");
stressRight = getState(domainRight, "stress");
stressRight(:, 1) = -1.0;
setState(domainRight, stressRight, "stress");

printTimes = 0:endTime;
output = OutState( ...
                  'outputFile', historyBase, ...
                  'printTimes', printTimes, ...
                  'matFileName', historyBase, ...
                  'saveHistory', 1, ...
                  'solvePrintTimes', 1);

solver = NonLinearImplicit( ...
                           'simulationparameters', simParam, ...
                           'domains', domains, ...
                           'interface', interfaces, ...
                           'output', output);
solver.simulationLoop();

% GReS stores P0 interface values on the slave grid (MortarSide.slave = 1).
centers = interfaces{1}.grids(1).surfaces.center;
end

function appendPlotStepMetadata(matFile, centers, direction, scenario, ...
                                method, augmentation, model)
if ~isfile(matFile)
    warning('Mortar:MissingHistory', ...
            'History MAT file was not found: %s', matFile);
    return
end
plotStepMeta = struct( ...
                      'model', model, ...
                      'centers', centers, ...
                      'verticalCoordinate', 3, ...
                      'transverseCoordinate', 2, ...
                      'direction', direction, ...
                      'scenario', scenario, ...
                      'method', method, ...
                      'augmentation', augmentation);
save(matFile, 'plotStepMeta', '-append');
end

function plotFixedComparison(records, config, direction)
fig = comparisonFigure(direction + " loading - fixed augmentation");
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, sprintf('%s loading, $c = 10^{%d}$', direction, ...
                      round(log10(config.fixedAugmentation))), ...
      'Interpreter', 'latex', 'FontSize', 18);

for iScenario = 1:numel(config.scenarios)
    ax = nexttile(layout, iScenario);
    values = nan(1, numel(config.methods));
    for iMethod = 1:numel(config.methods)
        rec = getRecord(records, direction, config.scenarios(iScenario), ...
                        config.methods(iMethod).id, config.fixedAugmentation);
        if ~isempty(rec) && rec.success
            values(iMethod) = rec.iterations;
        end
    end
    labels = string({config.methods.label});
    colors = vertcat(config.methods.color);
    drawBars(ax, values, labels, colors);
    title(ax, config.scenarioLabels(iScenario), 'Interpreter', 'latex', 'FontSize', 15);
    xlabel(ax, 'formulation', 'Interpreter', 'latex', 'FontSize', 14);
end

name = sprintf('%s_formulations_fixed_c.pdf', direction);
exportgraphics(fig, fullfile(config.figureDir, name), 'ContentType', 'vector');
end

function plotAugmentationSweep(records, config, direction, methods)

fig = comparisonFigure(direction + " loading - augmentation sweep");

layout = tiledlayout(fig, 1, 2, ...
                     'TileSpacing', 'compact', ...
                     'Padding', 'compact');

title(layout, ...
      direction + " loading - smooth-regularized and semi-smooth sweep", ...
      'Interpreter', 'latex', ...
      'FontSize', 18);

for iScenario = 1:numel(config.scenarios)
    ax = nexttile(layout, iScenario);

    % Rows correspond to augmentation values.
    % Columns correspond to the two formulations.
    values = nan(numel(config.augmentationSweep), numel(methods));

    for iMethod = 1:numel(methods)
        for iAugmentation = 1:numel(config.augmentationSweep)
            rec = getRecord(records, ...
                            direction, ...
                            config.scenarios(iScenario), ...
                            methods(iMethod).id, ...
                            config.augmentationSweep(iAugmentation));

            if ~isempty(rec) && rec.success
                values(iAugmentation, iMethod) = rec.iterations;
            end
        end
    end

    labels = compose('$10^{%d}$', ...
                     round(log10(config.augmentationSweep)));

    colors = vertcat(methods.color);
    methodLabels = string({methods.label});

    drawGroupedSweepBars(ax, ...
                         values, ...
                         labels, ...
                         colors, ...
                         methodLabels);

    title(ax, config.scenarioLabels(iScenario), ...
          'Interpreter', 'latex', ...
          'FontSize', 15);

    xlabel(ax, 'augmentation scaling $c$', ...
           'Interpreter', 'latex', ...
           'FontSize', 14);
end

name = sprintf('%s_augmentation_sweep.pdf', direction);

exportgraphics(fig, ...
               fullfile(config.figureDir, name), ...
               'ContentType', 'vector');

end

function drawGroupedSweepBars(ax, values, labels, colors, ...
                              methodLabels)

hold(ax, 'on');
box(ax, 'on');
grid(ax, 'on');

finiteValues = values(isfinite(values));

if isempty(finiteValues)
    yMax = 1;
else
    yMax = max(finiteValues);
end

dummyHeight = max(0.025 * yMax, 0.5);

nGroups = size(values, 1);
nMethods = size(values, 2);

% Narrow bars with a small space between the formulations.
barWidth = min(0.28, 0.64 / max(1, nMethods));

if nMethods == 1
    offsets = 0;
else
    offsets = ((1:nMethods) - (nMethods + 1) / 2) ...
      * (barWidth + 0.06);
end

barHandles = gobjects(1, nMethods);

for iMethod = 1:nMethods
    x = (1:nGroups) + offsets(iMethod);

    heights = values(:, iMethod).';
    failed = ~isfinite(heights);
    heights(failed) = dummyHeight;

    barHandles(iMethod) = bar(ax, ...
                              x, ...
                              heights, ...
                              barWidth, ...
                              'FaceColor', colors(iMethod, :), ...
                              'EdgeColor', [0.15, 0.15, 0.15], ...
                              'LineWidth', 0.8);

    for iGroup = 1:nGroups
        if failed(iGroup)
            plot(ax, ...
                 x(iGroup), ...
                 heights(iGroup) + 0.02 * yMax, ...
                 'rx', ...
                 'MarkerSize', 9, ...
                 'LineWidth', 1.8, ...
                 'HandleVisibility', 'off');

            text(ax, ...
                 x(iGroup), ...
                 heights(iGroup) + 0.055 * yMax, ...
                 'fail', ...
                 'Interpreter', 'latex', ...
                 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'bottom', ...
                 'FontSize', 10);

        end
    end
end

ax.XTick = 1:nGroups;
ax.XTickLabel = labels;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 11;

xlim(ax, [0.5, nGroups + 0.5]);
ylim(ax, [0, max(1, 1.22 * yMax)]);

ylabel(ax, 'total nonlinear iterations', ...
       'Interpreter', 'latex', ...
       'FontSize', 14);

legend(ax, ...
       barHandles, ...
       methodLabels, ...
       'Interpreter', 'latex', ...
       'Location', 'best', ...
       'FontSize', 11);

end

function rec = getRecord(records, direction, scenario, method, augmentation)
idx = findRecord(records, direction, scenario, method, augmentation);
if isempty(idx)
    rec = [];
else
    rec = records(idx);
end
end

function tag = scientificTag(value)
tag = replace(string(sprintf('%.0e', value)), ["+", "-"], ["p", "m"]);
end

function ensureDirectory(pathName)
if ~isfolder(pathName)
    mkdir(pathName);
end
end

function [bcLeft, bcRight] = setVerticalBC(Y, gridLeft, gridRight)
[bcLeft, bcRight] = baseBoundaryConditions(Y, gridLeft, gridRight);

bcRight.addBC('name', "z_load", 'type', "neumann", 'field', "surface", ...
              'variable', "displacements", 'entityListType', "tag", ...
              'entityList', 2, 'components', "z");
bcRight.addBCEvent("z_load", 'time', 0.0,  'value', 0.0);
bcRight.addBCEvent("z_load", 'time', 1.0,  'value', 0.0);
bcRight.addBCEvent("z_load", 'time', 6.0,  'value', -18.0);
bcRight.addBCEvent("z_load", 'time', 11.0, 'value', -18.0);
bcRight.addBCEvent("z_load", 'time', 16.0, 'value', 0.0);
end

function [bcLeft, bcRight] = setHorizontalBC(Y, gridLeft, gridRight)
[bcLeft, bcRight] = baseBoundaryConditions(Y, gridLeft, gridRight);

% This preserves the supplied horizontal benchmark: z traction is applied
% to surface tag 3 of the right block.
bcRight.addBC('name', "y_load", 'type', "neumann", 'field', "surface", ...
              'variable', "displacements", 'entityListType', "tag", ...
              'entityList', 3, 'components', "z");
bcRight.addBCEvent("y_load", 'time', 0.0,  'value', 0.0);
bcRight.addBCEvent("y_load", 'time', 1.0,  'value', 0.0);
bcRight.addBCEvent("y_load", 'time', 6.0,  'value', 5.0);
bcRight.addBCEvent("y_load", 'time', 11.0, 'value', 5.0);
bcRight.addBCEvent("y_load", 'time', 16.0, 'value', 0.0);
end

function [bcLeft, bcRight] = baseBoundaryConditions(Y, gridLeft, gridRight)
targetCoord = 0.5 * Y;
nLeft = find(abs(gridLeft.coordinates(:, 2) - targetCoord) < 1e-4 & ...
             abs(gridLeft.coordinates(:, 3)) < 1e-4);
nRight = find(abs(gridRight.coordinates(:, 2) - targetCoord) < 1e-4 & ...
              abs(gridRight.coordinates(:, 3)) < 1e-4);

bcLeft = Boundaries(gridLeft);
bcLeft.addBC('name', "fixBack", 'type', "dirichlet", 'field', "surface", ...
             'variable', "displacements", 'entityListType', "tag", ...
             'entityList', 5, 'components', "x");
bcLeft.addBCEvent("fixBack", 'time', 0.0, 'value', 0.0);
bcLeft.addBC('name', "y_bottom", 'type', "dirichlet", 'field', "node", ...
             'variable', "displacements", 'entityListType', "bcList", ...
             'entityList', nLeft, 'components', "y");
bcLeft.addBCEvent("y_bottom", 'time', 0.0, 'value', 0.0);
bcLeft.addBC('name', "z_bottom", 'type', "dirichlet", 'field', "surface", ...
             'variable', "displacements", 'entityListType', "tag", ...
             'entityList', 1, 'components', "z");
bcLeft.addBCEvent("z_bottom", 'time', 0.0, 'value', 0.0);

bcRight = Boundaries(gridRight);
bcRight.addBC('name', "x_load", 'type', "neumann", 'field', "surface", ...
              'variable', "displacements", 'entityListType', "tag", ...
              'entityList', 6, 'components', "x");
bcRight.addBCEvent("x_load", 'time', 0.0,  'value', 0.0);
bcRight.addBCEvent("x_load", 'time', 1.0,  'value', -5.0);
bcRight.addBCEvent("x_load", 'time', 6.0,  'value', -5.0);
bcRight.addBCEvent("x_load", 'time', 11.0, 'value', 0.0);
bcRight.addBCEvent("x_load", 'time', 16.0, 'value', 0.0);
bcRight.addBCEvent("x_load", 'time', 20.0, 'value', 1.0);
bcRight.addBC('name', "y_bottom", 'type', "dirichlet", 'field', "node", ...
              'variable', "displacements", 'entityListType', "bcList", ...
              'entityList', nRight, 'components', "y");
bcRight.addBCEvent("y_bottom", 'time', 0.0, 'value', 0.0);
bcRight.addBC('name', "z_bottom", 'type', "dirichlet", 'field', "surface", ...
              'variable', "displacements", 'entityListType', "tag", ...
              'entityList', 1, 'components', "z");
bcRight.addBCEvent("z_bottom", 'time', 0.0, 'value', 0.0);
end

function fig = comparisonFigure(figureName)

fig = figure( ...
             'Name', char(figureName), ...
             'Color', 'w', ...
             'Units', 'normalized', ...
             'Position', [0.1, 0.15, 0.8, 0.65]);

end

function drawBars(ax, values, labels, colors)
hold(ax, 'on');
box(ax, 'on');
grid(ax, 'on');

finiteValues = values(isfinite(values));
if isempty(finiteValues)
    yMax = 1;
else
    yMax = max(finiteValues);
end
dummyHeight = max(0.025 * yMax, 0.5);

barSpacing = 1.75;
xPositions = (1:numel(values)) * barSpacing;
barWidth = 0.55;

for i = 1:numel(values)
    failed = ~isfinite(values(i));
    height = values(i);
    if failed
        height = dummyHeight;
    end
    bar(ax, xPositions(i), height, 0.72, ...
        'FaceColor', colors(i, :), ...
        'EdgeColor', [0.15, 0.15, 0.15], ...
        'LineWidth', 0.8);

    if failed
        plot(ax, xPositions(i), height + 0.02 * yMax, 'rx', ...
             'MarkerSize', 9, 'LineWidth', 1.8);
        text(ax, xPositions(i), height + 0.055 * yMax, 'fail', ...
             'Interpreter', 'latex', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', 'FontSize', 11);
    end
end

ax.XTick = xPositions;
xlim(ax, [xPositions(1) - 0.7 * barSpacing, ...
          xPositions(end) + 0.7 * barSpacing]);
ax.XTickLabel = labels;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 11;
ylabel(ax, 'total nonlinear iterations', 'Interpreter', 'latex', 'FontSize', 14);
ylim(ax, [0, max(1, 1.22 * yMax)]);
end
