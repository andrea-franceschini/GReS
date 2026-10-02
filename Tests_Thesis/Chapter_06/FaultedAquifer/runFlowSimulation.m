function [grid, pressures, times] = runFlowSimulation(grid, meshInfo, config)
% Single-phase TPFA flow, on the native Grid. Time in years, pressure in kPa.
fprintf('\nFaulted aquifer: flow simulation\n');
materials = Materials('InputFlow/materials.xml');
bc = Boundaries(grid);
% Original far-field pressure conditions: x boundaries and free top.
outer = find(ismember(grid.surfaces.tag, [5 7 8]));
bc.addBC('name', "farField", 'type', "dirichlet", 'field', "surface", ...
    'variable', "pressure", 'entityListType', "bcList", 'entityList', outer);
bc.addBCEvent("farField", 'time', 0, 'value', 0);
bc.addBC('name', "wells", 'type', "dirichlet", 'field', "cell", ...
    'variable', "pressure", 'entityListType', "bcList", 'entityList', meshInfo.wellCells);
for k = 1:numel(config.wellTimes)
    bc.addBCEvent("wells", 'time', config.wellTimes(k), 'value', config.wellDrawdown(k));
end

domain = Discretizer('grid', grid, 'boundaries', bc, 'materials', materials);
domain.addPhysicsSolver('SinglePhaseFlowFVTPFA');
parameters = SimulationParameters('Start', 0, 'End', 10, 'DtInit', 1, 'DtMin', 1, ...
    'DtMax', 1, 'RelativeTolerance', 1e-8, 'AbsoluteTolerance', 1e-9, ...
    'MaxNLIteration', 10, 'LinearSolver', struct('useMatlab', 1));
output = OutState('outputFile', 'Outputs/outFlow', 'matFileName', 'Outputs/historyFlow', ...
    'printTimes', config.outputTimes, 'solvePrintTimes', 1);
solver = NonLinearImplicit('simulationparameters', parameters, 'domains', domain, 'output', output);
solver.simulationLoop();

% OutState now exposes results, not the obsolete matFile.pressure property.
times = [output.results.time];
pressures = [output.results.pressure];
% The current time loop may not print its initial state: supply the exact zero.
if isempty(times) || times(1) > 0
    times = [0 times];
    pressures = [zeros(grid.cells.num, 1) pressures];
end
assert(size(pressures, 1) == grid.cells.num, 'Pressure/cell ordering mismatch.');
assert(all(isfinite(pressures), 'all'), 'Nonfinite flow solution.');
assert(abs(times(end) - 10) < 1e-10, 'Flow simulation did not reach 10 years.');
save('Outputs/pressureTransfer.mat', 'pressures', 'times', 'meshInfo');
wellTable = array2table([times(:) pressures(meshInfo.wellCells, :)'], ...
    'VariableNames', {'time_years', 'well1_drawdown_kPa', 'well2_drawdown_kPa'});
writetable(wellTable, 'Outputs/wellDrawdown.csv');
end
