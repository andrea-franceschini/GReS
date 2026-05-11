function [maxF, avgP, avgQ, ratio] = simulate(solver, forceX, forceZ, druckerPrager)

    % scale boundary conditions
    solver.domains.bcs.scaleBC('load_top',forceZ);
    solver.domains.bcs.scaleBC('load_lateral',forceX);

    solver.simulationLoop();
    state = solver.domains.getState();
    stress = state.stress;

    [maxF, avgP, avgQ, ratio] = ...
        check_drucker_prager(solver.domains.grid, stress, druckerPrager);

end
