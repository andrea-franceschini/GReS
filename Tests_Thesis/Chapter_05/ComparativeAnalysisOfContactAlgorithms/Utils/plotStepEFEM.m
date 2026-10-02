function [fig, profile] = plotStep(matFile, timeStep)
% PLOTSTEP Plot EFEM contact profiles along the crack central vertical axis.
%   plotStep(MATFILE, TIMESTEP) reads a GReS history MAT file and plots
%   t_N, ||t_T||_2 and ||g_T||_2 at the requested saved time/load step.
%   The figure is exported beside MATFILE. TIMESTEP is matched against the
%   saved output.time values first; if no time matches, it is used as a
%   one-based history index.

arguments
    matFile (1, 1) string
    timeStep (1, 1) double {mustBeFinite, mustBeNonnegative}
end

if ~isfile(matFile)
    error('plotStep:MissingFile', 'MAT file not found: %s', matFile);
end

history = load(matFile);
if ~isfield(history, 'output') || isempty(history.output)
    error('plotStep:MissingOutput', ...
      'The MAT file does not contain the GReS variable "output".');
end

[step, stepIndex, savedTime] = selectStep(history.output, timeStep);
[tN, tTNorm, gTNorm] = extractEFEMFields(step);
nValues = numel(tN);

if isfield(history, 'plotStepMeta') && ...
    isfield(history.plotStepMeta, 'centers') && ...
    ~isempty(history.plotStepMeta.centers)
    centers = history.plotStepMeta.centers;
else
    centers = inferCenters(nValues);
    warning('plotStep:InferredGeometry', ...
      ['plotStepMeta was not present. Centers were reconstructed for the ', ...
      'supplied 10-by-10 crack grid.']);
end

if size(centers, 1) ~= nValues
    error('plotStep:GeometrySize', ...
      'The history has %d crack values but %d crack centers.', ...
      nValues, size(centers, 1));
end

[z, tN, tTNorm, gTNorm] = centralAxisProfile( ...
  centers, tN, tTNorm, gTNorm);

profile = table(z, tN, tTNorm, gTNorm, ...
  'VariableNames', {'z', 'tN', 'tTNorm', 'gTNorm'});

fig = figure( ...
  'Name', sprintf('EFEM step %d', stepIndex), ...
  'Color', 'w', ...
  'Units', 'pixels', ...
  'Position', [100, 100, 1300, 520]);
layout = tiledlayout(fig, 1, 3, ...
  'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, sprintf('EFEM - saved time $t = %g$', savedTime), ...
  'Interpreter', 'latex', 'FontSize', 17);

plotProfile(nexttile(layout, 1), tN, z, '$t_N$');
plotProfile(nexttile(layout, 2), tTNorm, z, '$\|\mathbf{t}_T\|_2$');
plotProfile(nexttile(layout, 3), gTNorm, z, '$\|\mathbf{g}_T\|_2$');

[folder, baseName] = fileparts(matFile);
outFile = fullfile(folder, sprintf('%s_step_%03d_profiles.png', ...
  baseName, stepIndex));
exportgraphics(fig, outFile, 'Resolution', 300);
fprintf('Profile figure saved to: %s\n', outFile);
end

function [step, index, savedTime] = selectStep(output, requested)
times = nan(1, numel(output));
for i = 1:numel(output)
    if isfield(output(i), 'time') && isscalar(output(i).time)
        times(i) = output(i).time;
    end
end

tolerance = 100 * eps(max(1, abs(requested)));
index = find(abs(times - requested) <= tolerance, 1, 'first');
if isempty(index)
    if requested == round(requested) && requested >= 1 && requested <= numel(output)
        index = round(requested);
    else
        available = strjoin(compose('%g', times(isfinite(times))), ', ');
        error('plotStep:UnknownStep', ...
          'No saved time %g. Available saved times: %s', requested, available);
    end
end

step = output(index);
if isfield(step, 'time') && isscalar(step.time)
    savedTime = step.time;
else
    savedTime = index;
end
end

function [tN, tTNorm, gTNorm] = extractEFEMFields(step)
if ~isfield(step, 'traction')
    error('plotStep:MissingTraction', ...
      'The selected EFEM output does not contain "traction".');
end
if ~isfield(step, 'fractureJump')
    error('plotStep:MissingJump', ...
      'The selected EFEM output does not contain "fractureJump".');
end

traction = step.traction(:);
jump = step.fractureJump(:);
if mod(numel(traction), 3) ~= 0 || mod(numel(jump), 3) ~= 0
    error('plotStep:FieldShape', ...
      'EFEM traction and fractureJump must contain interleaved 3-vectors.');
end

traction = reshape(traction, 3, []);
jump = reshape(jump, 3, []);
tN = traction(1, :).';
tTNorm = vecnorm(traction(2:3, :), 2, 1).';
gTNorm = vecnorm(jump(2:3, :), 2, 1).';
end

function centers = inferCenters(nValues)
nSide = round(sqrt(nValues));
if nSide^2 ~= nValues
    error('plotStep:MissingGeometry', ...
      ['This older MAT file has no plotStepMeta and its %d values cannot ', ...
      'be mapped to the supplied square crack grid.'], nValues);
end
y = ((1:nSide) - 0.5) * (10 / nSide);
z = ((1:nSide) - 0.5) * (15 / nSide);
[Y, Z] = ndgrid(y, z);
centers = [2.5 * ones(nValues, 1), Y(:), Z(:)];
end

function [z, tN, tTNorm, gTNorm] = centralAxisProfile( ...
  centers, tN, tTNorm, gTNorm)

y = centers(:, 2);
zAll = centers(:, 3);
yMid = 0.5 * (min(y) + max(y));
distance = abs(y - yMid);
minDistance = min(distance);
tolerance = max(1e-10, 1e-8 * max(1, range(y)));
central = distance <= minDistance + tolerance;

zSelected = zAll(central);
tNSelected = tN(central);
tTSelected = tTNorm(central);
gTSelected = gTNorm(central);

roundedZ = round(zSelected, 10);
[z, ~, group] = unique(roundedZ, 'sorted');
tN = accumarray(group, tNSelected, [], @mean);
tTNorm = accumarray(group, tTSelected, [], @mean);
gTNorm = accumarray(group, gTSelected, [], @mean);
end

function plotProfile(ax, values, z, xLabel)
plot(ax, values, z, '-o', ...
  'Color', [0.850, 0.125, 0.098], ...
  'MarkerFaceColor', [0.850, 0.125, 0.098], ...
  'MarkerEdgeColor', [0.850, 0.125, 0.098], ...
  'LineWidth', 1.8, ...
  'MarkerSize', 5);
box(ax, 'on');
grid(ax, 'on');
xlabel(ax, xLabel, 'Interpreter', 'latex', 'FontSize', 15);
ylabel(ax, '$z$', 'Interpreter', 'latex', 'FontSize', 15);
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 12;
ylim(ax, [min(z), max(z)]);
end
