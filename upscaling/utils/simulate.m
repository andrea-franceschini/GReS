function [maxF, avgP, avgQ, ratio] = simulate(fea, forceX, forceZ, druckerPrager)

    % scale boundary conditions
    fea.solver.domains.bcs.scaleBC('load_top',forceZ);
    fea.solver.domains.bcs.scaleBC('load_lateral',forceX);

    % reset the state
    fea.solver.domains.state = copy(fea.initState);
    fea.solver.domains.stateOld = copy(fea.initState);

    fea.solver.simulationLoop();
    stress = fea.solver.domains.state.data.stress;

    [maxF, avgP, avgQ, ratio] = ...
        check_drucker_prager(fea.solver.domains.grid, stress, druckerPrager);

end
