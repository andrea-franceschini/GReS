function [avgF3, avgP, avgQ] = solvePQ(solver, forceX, druckerPrager, tol)

    forceZ1 = forceX;
    [avgF1, ~, ~, ratio1] = simulate(solver, forceX, forceZ1, druckerPrager);
    fprintf(' -- %5d%15.6e%15.6e%15.6e\n', 0, forceZ1, avgF1, ratio1);

    forceZ2 = 4.0 * forceX;
    [avgF2, ~, ~, ratio2] = simulate(solver, forceX, forceZ2, druckerPrager);
    fprintf(' -- %5d%15.6e%15.6e%15.6e\n', 1, forceZ2, avgF2, ratio2);

    if avgF1 > 1.0 || avgF2 < 1.0
        error('Bad setting!');
    end

    avgF3 = 0.0;
    ratio3 = 0.0;
    it = 1;

    while (forceZ2 - forceZ1 > tol * forceX) && abs(ratio3 - 0.5) > tol

        forceZ3 = (forceZ1 + forceZ2) / 2.0;
        [avgF3, avgP, avgQ, ratio3] = simulate(solver, forceX, forceZ3, druckerPrager);

        if (ratio3 - 0.5) * (ratio1 - 0.5) > 0
            forceZ1 = forceZ3;
        else
            forceZ2 = forceZ3;
        end

        it = it + 1;
        fprintf(' -- %5d%15.6e%15.6e%15.6e\n', it, forceZ3, avgF3, ratio3);

    end

end
