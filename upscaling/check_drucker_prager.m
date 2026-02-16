function [maxF, avgP, avgQ, ratio] = ...
    check_drucker_prager(grid, stress, druckerPrager)

    % Volumes
    vols = grid.topology.cellVolume;

    % Regions
    regions = grid.topology.cellTag;

    % Zone IDs
    ID_DZ1  = regions == druckerPrager.DamageZone1.zoneID;
    ID_core = regions == druckerPrager.Core.zoneID;
    ID_DZ2  = regions == druckerPrager.DamageZone2.zoneID;

    % Mean properties
    cohes_DZ1  = druckerPrager.DamageZone1.cohesion_mean;
    phi_DZ1    = druckerPrager.DamageZone1.friction_angle_mean;

    cohes_core = druckerPrager.Core.cohesion_mean;
    phi_core   = druckerPrager.Core.friction_angle_mean;

    cohes_DZ2  = druckerPrager.DamageZone2.cohesion_mean;
    phi_DZ2    = druckerPrager.DamageZone2.friction_angle_mean;

    % Randomize cohesion and friction angle
    phi_DZ1  = generate_truncated_random_numbers(phi_DZ1, ...
                druckerPrager.DamageZone1.friction_angle_std, 0, inf, ...
                sum(ID_DZ1));
    phi_core = generate_truncated_random_numbers(phi_core, ...
                druckerPrager.Core.friction_angle_std, 0, inf, ...
                sum(ID_core));
    phi_DZ2  = generate_truncated_random_numbers(phi_DZ2, ...
                druckerPrager.DamageZone2.friction_angle_std, 0, inf, ...
                sum(ID_DZ2));

    cohes_DZ1  = generate_truncated_random_numbers(cohes_DZ1, ...
                   druckerPrager.DamageZone1.cohesion_std, 0, inf, ...
                   sum(ID_DZ1));
    cohes_core = generate_truncated_random_numbers(cohes_core, ...
                   druckerPrager.Core.cohesion_std, 0, inf, ...
                   sum(ID_core));
    cohes_DZ2  = generate_truncated_random_numbers(cohes_DZ2, ...
                   druckerPrager.DamageZone2.cohesion_std, 0, inf, ...
                   sum(ID_DZ2));

    % ---- Convert to Drucker-Prager parameters ----

    phi_DZ1  = deg2rad(phi_DZ1);
    phi_core = deg2rad(phi_core);
    phi_DZ2  = deg2rad(phi_DZ2);

    A_DZ1 = 3*tan(phi_DZ1) ./ sqrt(9 + 12*tan(phi_DZ1).^2);
    B_DZ1 = 3*cohes_DZ1 ./ sqrt(9 + 12*tan(phi_DZ1).^2);

    A_core = 3*tan(phi_core) ./ sqrt(9 + 12*tan(phi_core).^2);
    B_core = 3*cohes_core ./ sqrt(9 + 12*tan(phi_core).^2);

    A_DZ2 = 3*tan(phi_DZ2) ./ sqrt(9 + 12*tan(phi_DZ2).^2);
    B_DZ2 = 3*cohes_DZ2 ./ sqrt(9 + 12*tan(phi_DZ2).^2);

    % ---- Stress components (Voigt: xx, yy, zz, yz, xz, xy) ----

    sigma_x = stress(:,1);
    sigma_y = stress(:,2);
    sigma_z = stress(:,3);
    tau_yz  = stress(:,4);
    tau_xz  = stress(:,5);
    tau_xy  = stress(:,6);

    % ---- Mean pressure and deviatoric stress ----

    p = (sigma_x + sigma_y + sigma_z) / 3.0;

    q = sqrt( ...
        ((sigma_x - p).^2 + ...
         (sigma_y - p).^2 + ...
         (sigma_z - p).^2)/2.0 + ...
         tau_xy.^2 + tau_yz.^2 + tau_xz.^2 );

    % ---- Yield function ----

    f_DZ1  = q(ID_DZ1)  ./ (-A_DZ1  .* p(ID_DZ1)  + B_DZ1);
    f_core = q(ID_core) ./ (-A_core .* p(ID_core) + B_core);
    f_DZ2  = q(ID_DZ2)  ./ (-A_DZ2  .* p(ID_DZ2)  + B_DZ2);

    F_DZ1  = sum(f_DZ1  .* vols(ID_DZ1))  / sum(vols(ID_DZ1));
    F_core = sum(f_core .* vols(ID_core)) / sum(vols(ID_core));
    F_DZ2  = sum(f_DZ2  .* vols(ID_DZ2))  / sum(vols(ID_DZ2));

    % ---- Plastic volume ratio ----

    vols_DZ1 = vols(ID_DZ1);
    vols_core = vols(ID_core);
    vols_DZ2 = vols(ID_DZ2);
    plastic_vol = sum(vols_DZ1(f_DZ1 > 1)) + ...
        sum(vols_core(f_core > 1)) + sum(vols_DZ2(f_DZ2 > 1));

    total_vol = sum(vols(ID_DZ1)) + sum(vols(ID_core)) + sum(vols(ID_DZ2));

    ratio = plastic_vol / total_vol;

    % ---- Average p and q ----

    avgP = ( sum(p(ID_DZ1).*vols_DZ1) + ...
             sum(p(ID_core).*vols_core) + ...
             sum(p(ID_DZ2).*vols_DZ2) ) / total_vol;

    avgQ = ( sum(q(ID_DZ1).*vols_DZ1) + ...
             sum(q(ID_core).*vols_core) + ...
             sum(q(ID_DZ2).*vols_DZ2) ) / total_vol;

    maxF = max([F_DZ1, F_core, F_DZ2]);

end
