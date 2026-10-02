function result = runContactMechanicsSimulation(grid, drawdown, times, config)

% One-way mechanics driven by the flow history on exactly the same cells.
% E and stresses in MPa. Drawdown is positive; pore-pressure change is negative.
fprintf('\nFaulted aquifer: contact mechanics simulation\n');
assert(size(drawdown, 1) == grid.cells.num && size(drawdown, 2) == numel(times), ...
    'Flow history must contain one row per mechanical cell.');
materials = Materials('InputMech/materials.xml');

bc = Boundaries(grid);
% Roller constraints: zero normal displacement on sides and bottom;
% the top remains traction-free, including at the fault intersection.
faceTags = {3, [4 6], [5 7]};
components = [3 2 1];
for k = 1:3
    nodes = unique(grid.surfaces.connectivity(ismember(grid.surfaces.tag, faceTags{k}), :));
    name = "roller" + k;
    bc.addBC('name', name, 'type', "dirichlet", 'field', "node", ...
        'variable', "displacements", 'entityListType', "bcList", ...
        'entityList', nodes, 'components', components(k));
    bc.addBCEvent(name, 'time', 0, 'value', 0);
end

% Consistent pressure-to-force map, integrated with native GReS elements.
% This is the same weak load Q*delta_p as Poromechanics' volume-force BC.
% Register its components as nodal Neumann data, without generating BC files.
Q = pressureLoadMatrix(grid);
forces = Q * (-1e-3 * drawdown);       % kPa drawdown -> negative MPa pressure change
nodes = (1:grid.nNodes)';

for component = 1:3
    name = "pressureLoad" + component;
    bc.addBC('name', name, 'type', "neumann", 'field', "node", ...
        'variable', "displacements", 'entityListType', "bcList", ...
        'entityList', nodes, 'components', component);
    for k = 1:numel(times)
        bc.addBCEvent(name, 'time', times(k), 'value', forces(component:3:end, k));
    end
end

domain = Discretizer('grid', grid, 'boundaries', bc, 'materials', materials);
domain.addPhysicsSolver('Poromechanics');
initializeFaultStress(domain, config);

% Master is stiff rock (tag 1); slave is the layered soil (tag 2).
contactInput = struct('masterDomain', 1, 'slaveDomain', 1, ...
    'masterSurface', 1, 'slaveSurface', 2, 'multiplierType', "P0", ...
    'stabilizationScale', config.stabilizationScale, ...
    'Quadrature', struct('type', "SegmentBasedQuadrature", 'gaussOrder', 4), ...
    'Coulomb', struct('cohesion', config.cohesion, 'frictionAngle', config.frictionAngle), ...
    'ActiveSet', struct('resetActiveSet', 0, 'forceStickBoundary', "z", ...
        'Tolerances', struct('sliding', 1e-6, 'normalGap', 1e-2, ...
        'normalTraction', 1e-6, 'tangentialViolation', 1e-3, ...
        'minLimitTraction', 0, 'areaChange', 1e-4, 'maxStateChange', 20)));
interfaces = InterfaceSolver.add('SolidMechanicsContactNew', domain, contactInput);
% GReS initializes fault traction from the initial bulk stress in its local
% frame, and pins slave faces touching the bottom displacement boundary.
parameters = SimulationParameters('Start', 0, 'End', 10, 'DtInit', 1, 'DtMin', 1e-4, ...
    'DtMax', 1, 'incrementFactor', 1.25, 'choppingFactor', 2, ...
    'RelativeTolerance', 1e-9, 'AbsoluteTolerance', 1e-8, 'MaxNLIteration', 10, ...
    'MaxConfigurationIteration', 20, 'LinearSolver', struct('useMatlab', 1));

output = OutState('outputFile', 'Outputs/outMech', 'matFileName', 'Outputs/historyMech', ...
    'printTimes', config.outputTimes, 'solvePrintTimes', 1);

solver = NonLinearImplicit('simulationparameters', parameters, 'domains', domain, ...
    'interface', interfaces, 'output', output);

solver.simulationLoop();

result.displacements = domain.getState('displacements');
result.contact = interfaces{1}.getState();
result.activeSet = interfaces{1}.activeSet.curr;
assert(all(isfinite(result.displacements)), 'Nonfinite mechanical solution.');
assert(all(isfinite(result.contact.traction)), 'Nonfinite fault traction.');
result.maxDisplacement = max(vecnorm(reshape(result.displacements, 3, [])', 2, 2));
save('Outputs/mechanicalResult.mat', 'result', 'times');
end
