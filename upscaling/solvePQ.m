function [avgF, avgP, avgQ] = solvePQ(fea, forceX, druckerPrager, tol)

    [forceZ1, forceZ2] = initial_guess(fea, forceX, druckerPrager);

    target = 0.5;

    % --- Initial bracket ---
    a = forceZ1;
    [~, ~, ~, ratio_a] = simulate(fea, forceX, a, druckerPrager);
    fa = ratio_a - target;

    b = forceZ2;
    [~, ~, ~, ratio_b] = simulate(fea, forceX, b, druckerPrager);
    fb = ratio_b - target;

    if fa * fb > 0
        error('Bad initial bracket!');
    end

    % --- Ensure ordering a < b ---
    if a > b
        [a,b] = deal(b,a);
        [fa,fb] = deal(fb,fa);
    end

    % Make |fb| <= |fa|
    if abs(fa) < abs(fb)
        [a,b] = deal(b,a);
        [fa,fb] = deal(fb,fa);
    end

    c = a;
    fc = fa;

    d = b - a;
    e = d;

    it = 0;

    while abs(b - a) > tol * forceX && abs(fb) > tol

        if fa ~= fc && fb ~= fc
            % Inverse quadratic interpolation
            s = a*fb*fc/((fa-fb)*(fa-fc)) ...
              + b*fa*fc/((fb-fa)*(fb-fc)) ...
              + c*fa*fb/((fc-fa)*(fc-fb));
        else
            % Secant step
            s = b - fb*(b-a)/(fb-fa);
        end

        % --- Acceptance conditions ---
        cond1 = (s <= min(a,b)) || (s >= max(a,b));
        cond2 = abs(e) < tol;
        cond3 = abs(s-b) >= abs(b-c)/2;

        if cond1 || cond2 || cond3
            % Bisection fallback
            s = (a + b)/2;
            e = b - a;
        else
            e = d;
        end

        d = b - s;

        [avgF, avgP, avgQ, ratio_s] = simulate(fea, forceX, s, druckerPrager);
        fs = ratio_s - target;

        c = a;
        fc = fa;

        if fa * fs < 0
            b = s;
            fb = fs;
        else
            a = s;
            fa = fs;
        end

        % --- Ensure ordering again ---
        if a > b
            [a,b] = deal(b,a);
            [fa,fb] = deal(fb,fa);
        end

        % Make |fb| <= |fa|
        if abs(fa) < abs(fb)
            [a,b] = deal(b,a);
            [fa,fb] = deal(fb,fa);
        end

        it = it + 1;
        fprintf(' -- %5d%15.6e%15.6e%15.6e\n', it, s, avgF, ratio_s);
    end

end

function [forceZ1, forceZ2] = initial_guess(fea, forceX, druckerPrager)
    % Volumes
    vols = fea.solver.domains.grid.topology.cellVolume;

    % Regions
    regions = fea.solver.domains.grid.topology.cellTag;

    % Zone IDs
    ID_DZ1  = regions == druckerPrager.DamageZone1.zoneID;
    ID_core = regions == druckerPrager.Core.zoneID;
    ID_DZ2  = regions == druckerPrager.DamageZone2.zoneID;

    vols_DZ1 = sum(vols(ID_DZ1));
    vols_core = sum(vols(ID_core));
    vols_DZ2 = sum(vols(ID_DZ2));

    cohes_DZ1  = druckerPrager.DamageZone1.cohesion_mean;
    phi_DZ1    = druckerPrager.DamageZone1.friction_angle_mean;

    cohes_core = druckerPrager.Core.cohesion_mean;
    phi_core   = druckerPrager.Core.friction_angle_mean;

    cohes_DZ2  = druckerPrager.DamageZone2.cohesion_mean;
    phi_DZ2    = druckerPrager.DamageZone2.friction_angle_mean;

    cohes = cohes_DZ1*vols_DZ1 + cohes_core*vols_core + cohes_DZ2*vols_DZ2;
    phi = phi_DZ1*vols_DZ1 + phi_core*vols_core + phi_DZ2*vols_DZ2;
    cohes = cohes / (vols_DZ1 + vols_core + vols_DZ2);
    phi = deg2rad(phi) / (vols_DZ1 + vols_core + vols_DZ2);
    A = 3*tan(phi) / sqrt(9 + 12*tan(phi)^2);
    B = 3*cohes / sqrt(9 + 12*tan(phi)^2);

    % Since p = (2*fx+fz)/3 and q = (fz-fx)/sqrt(3)
    % we have fx = p-q/sqrt(3) and fz = p+2*q/sqrt(3)
    % moreover q = A*p+B (with averaged values), thus:
    % fx = p*(1-A/sqrt(3))-B/sqrt(3) and p = (fx+B/sqrt(3))/(1-A/sqrt(3))
    % Finally fz = p*(1+2*A/sqrt(3))+2*B/sqrt(3)
    forceZ1 = 0.8*(((1+2*A/sqrt(3))*(forceX + B/sqrt(3)))/(1 - A/sqrt(3)) + 2*B/sqrt(3));
    forceZ2 = 1.2*(((1+2*A/sqrt(3))*(forceX + B/sqrt(3)))/(1 - A/sqrt(3)) + 2*B/sqrt(3));
end