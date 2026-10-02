function result = main(varargin)
% Poisson problems: curved interface (thesis section 4.4.2).
if nargin == 1 && isstruct(varargin{1}) && isfield(varargin{1}, 'SuiteRoot')
    result = runExperiment(varargin{1});
    return
end
suiteRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(suiteRoot, 'Utils'));
result = thesisRun(suiteRoot, 'Cases', "PoissonProblemsCurvedInterface", varargin{:});
end

function result = runExperiment(options)
% Poisson curved interface: independent mesh, boundary conditions, and solver setup.
cfg = jsondecode(fileread('Input/config.json'));
if options.Smoke
    cfg.leftSizes = cfg.leftSizes(1);
    cfg.rightSizes = cfg.rightSizes(1);
    cfg.orders = 1;
    cfg.quadratureOrders = cfg.quadratureOrders(1);
end

% Exact solution and source term.
u = @(x, y, z) cos(pi * y) .* cos(pi * z) .* (2 * x - x.^2 + sin(pi * x));
f = @(x, y, z) cos(pi * y) .* cos(pi * z) .* ...
    (-2 - 3 * pi^2 * sin(pi * x) - 4 * pi^2 * x + 2 * pi^2 * x.^2);
ux = @(x, y, z) cos(pi * y) .* cos(pi * z) .* (2 - 2 * x + pi * cos(pi * x));
uy = @(x, y, z) -pi * sin(pi * y) .* cos(pi * z) .* (2 * x - x.^2 + sin(pi * x));
uz = @(x, y, z) -pi * cos(pi * y) .* sin(pi * z) .* (2 * x - x.^2 + sin(pi * x));
methods = ["SegmentBasedQuadrature", "ElementBasedQuadrature", "RBFquadrature"];
rows = [];

for order = reshape(cfg.orders, 1, [])
    for method = methods
        qOrders = reshape(cfg.quadratureOrders, 1, []);
        if method == "SegmentBasedQuadrature"
            qOrders = 5;
        end
        for qOrder = qOrders
            for level = 1:numel(cfg.leftSizes)
                left = cfg.leftSizes(level);
                right = cfg.rightSizes(level);
                domains = Discretizer.empty;
                domains(1) = makeDomain(1, left, order, cfg, u, f, ux, uy, uz, [1, 2, 3, 4, 5]); %#ok<AGROW>
                domains(2) = makeDomain(2, right, order, cfg, u, f, ux, uy, uz, [1, 2, 3, 4, 6]); %#ok<AGROW>

                interfaceFaces = [6 5];
                quadrature = struct('type', method, 'gaussOrder', qOrder, 'nInt', cfg.nInt);
                interfaceInput = struct('masterDomain', 1, 'slaveDomain', 2, ...
                    'masterSurface', interfaceFaces(1), 'slaveSurface', interfaceFaces(2), ...
                    'multiplierType', string(cfg.multiplierType), 'Quadrature', quadrature);
                interfaces = InterfaceSolver.add('MeshTying', domains, interfaceInput);

                label = sprintf('%s_Q%d_order%d_level%d', method, order, qOrder, level);
                simulation = SimulationParameters('Start', 0, 'End', 1, 'DtInit', 1, ...
                    'DtMin', 1, 'DtMax', 1, 'RelativeTolerance', 1e-10, ...
                    'AbsoluteTolerance', 1e-11, 'MaxNLIteration', 20);
                output = OutState('outputFile', fullfile('Output', label), 'printTimes', 1);
                solver = NonLinearImplicit('simulationparameters', simulation, ...
                    'domains', domains, 'interface', interfaces, 'output', output);
                solver.simulationLoop();

                l2Squared = 0;
                h1Squared = 0;
                maxNodalError = 0;
                for domainIndex = 1:2
                    poisson = domains(domainIndex).getPhysicsSolver('Poisson');
                    [l2, h1] = poisson.computeError();
                    l2Squared = l2Squared + l2^2;
                    h1Squared = h1Squared + h1^2;
                    maxNodalError = max(maxNodalError, max(abs(poisson.getState('err'))));
                end
                row = struct('method', method, 'elementOrder', order, ...
                    'quadratureOrder', qOrder, 'level', level, 'h', 1 / max(left, right), ...
                    'L2', sqrt(l2Squared), 'H1', sqrt(h1Squared), ...
                    'maxNodalError', maxNodalError, 'maxMultiplierError', NaN);
                rows = [rows; row]; %#ok<AGROW>
            end
        end
    end
end

result.errors = struct2table(rows);
result.errors.L2rate = nan(height(result.errors), 1);
result.errors.H1rate = nan(height(result.errors), 1);
for i = 2:height(result.errors)
    previous = rows(i - 1);
    current = rows(i);
    sameStudy = previous.method == current.method && ...
        previous.elementOrder == current.elementOrder && ...
        previous.quadratureOrder == current.quadratureOrder;
    if ~sameStudy || current.h >= previous.h
        continue
    end
    ratio = log(previous.h / current.h);
    result.errors.L2rate(i) = log(previous.L2 / current.L2) / ratio;
    result.errors.H1rate(i) = log(previous.H1 / current.H1) / ratio;
end
writetable(result.errors, fullfile('Output', 'errors.csv'));
plotConvergence(result.errors, cfg, methods);
end

function domain = makeDomain(index, n, order, cfg, u, f, ux, uy, uz, fixedFaces)
% Geometry and exterior Dirichlet faces for this experiment.
bounds = [index - 1 index; 0 1; 0 1];

grid = thesisHexMesh([n n n], bounds, order);
% Map both sides to the same curved x=1 interface.
coords = grid.coordinates;
shift = cfg.curvatureAmplitude * sin(pi * coords(:, 2)) .* sin(pi * coords(:, 3));

if index == 1
    grid.coordinates(:, 1) = coords(:, 1) .* (1 + shift);
else
    grid.coordinates(:, 1) = 1 + shift + (coords(:, 1) - 1) .* (1 - shift);
end

boundary = Boundaries(grid);
faces = ismember(grid.surfaces.tag, fixedFaces);
nodes = unique(grid.surfaces.connectivity(faces, :));
coords = grid.coordinates(nodes, :);
values = u(coords(:, 1), coords(:, 2), coords(:, 3));
thesisAddBC(boundary, "exact", "dirichlet", "node", "u", nodes, 1, values);

domain = Discretizer('grid', grid, 'boundaries', boundary);
domain.addPhysicsSolver('Poisson', 'gaussOrder', max(4, 2 * order + 2));
domain.getPhysicsSolver('Poisson').setAnalSolution(u, f, ux, uy, uz);
end

function plotConvergence(errors, cfg, methods)

fig = figure('Color', 'w');
fig.Position(3:4) = [1000 520];

tl = tiledlayout(1, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

lineStyles = {'-', '--', ':', '-.'};
markers    = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h'};

lineWidth  = 1.8;
markerSize = 6;

quantities = ["L2", "H1"];
yLabels = { ...
    '$\|e\|_{L^2}$', ...
    '$\|e\|_{H^1}$'};

% Store handles only once for the common legend
legendHandles = gobjects(0);
legendLabels  = {};

for k = 1:2

    ax = nexttile;
    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'on');

    quantity = quantities(k);

    for order = reshape(cfg.orders, 1, [])

        for method = methods

            qOrders = unique( ...
                errors.quadratureOrder(errors.method == method));

            for iq = 1:numel(qOrders)

                q = qOrders(iq);

                selected = ...
                    errors.elementOrder == order & ...
                    errors.method == method & ...
                    errors.quadratureOrder == q;

                if ~any(selected)
                    continue
                end

                h   = errors.h(selected);
                err = errors.(quantity)(selected);

                [h, idx] = sort(h);
                err = err(idx);

                p = loglog(ax, h, err, ...
                    'LineStyle', ...
                        lineStyles{mod(iq - 1, numel(lineStyles)) + 1}, ...
                    'Marker', ...
                        markers{mod(order - 1, numel(markers)) + 1}, ...
                    'LineWidth', lineWidth, ...
                    'MarkerSize', markerSize);

                % Collect legend entries only from first panel
                if k == 1
                    legendHandles(end + 1) = p;
                    legendLabels{end + 1} = sprintf( ...
                        '%s, $q=%d$, $p=%d$', ...
                        method, q, order);
                end
            end
        end
    end

    xlabel(ax, '$h$', 'Interpreter', 'latex');
    ylabel(ax, yLabels{k}, 'Interpreter', 'latex');

    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 12;
    ax.LineWidth = 0.8;
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';

end

% One common legend BELOW the plots
lgd = legend(legendHandles, legendLabels, ...
    'Interpreter', 'latex', ...
    'Orientation', 'horizontal', ...
    'NumColumns', 3, ...
    'FontSize', 10, ...
    'Box', 'off');

lgd.Layout.Tile = 'south';

exportgraphics(fig, ...
    fullfile('Output', 'convergence.pdf'), ...
    'ContentType', 'vector');

end
