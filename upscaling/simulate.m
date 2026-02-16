function [maxF, avgP, avgQ, ratio] = simulate(fea, forceX, forceZ, druckerPrager)

    axial_load = fea.solver.domains.bcs.db('load_top');
    values = axial_load.data.availVals(:,1);
    axial_load.data.availVals(:,1) = values * forceZ;
    fea.solver.domains.bcs.db('load_top') = axial_load;

    confining_load = fea.solver.domains.bcs.db('load_lateral');
    values = confining_load.data.availVals(:,1);
    confining_load.data.availVals(:,1) = values * forceX;
    fea.solver.domains.bcs.db('load_lateral') = confining_load;

    % reset the state
    fea.solver.domains.state = copy(fea.initState);
    fea.solver.domains.stateOld = copy(fea.initState);

    fea.solver.simulationLoop();
    stress = fea.solver.domains.state.data.stress;

    axial_load = fea.solver.domains.bcs.db('load_top');
    values = axial_load.data.availVals(:,1);
    axial_load.data.availVals(:,1) = values / forceZ;
    fea.solver.domains.bcs.db('load_top') = axial_load;

    confining_load = fea.solver.domains.bcs.db('load_lateral');
    values = confining_load.data.availVals(:,1);
    confining_load.data.availVals(:,1) = values / forceX;
    fea.solver.domains.bcs.db('load_lateral') = confining_load;

    [maxF, avgP, avgQ, ratio] = ...
        check_drucker_prager(fea.solver.domains.grid, stress, druckerPrager);

end
