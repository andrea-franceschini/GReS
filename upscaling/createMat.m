function stiffnesses = createMat(grid, meshProp, rock, blockSize, jointFamilies)

    % Get joint/rock parameters
    [E_rock, nu_rock, azimuths, dips, Kn_values, Ks_values, ...
        numbers_of_joints] = set_parameters(rock, jointFamilies);

    H = meshProp.H;
    theta = deg2rad(meshProp.fault_angle);
    cos_t = cos(theta);
    sin_t = sin(theta);

    pt   = [0, 0, -H/2];
    norm = [-sin_t, 0, cos_t];
    d    = -dot(norm, pt);

    coords = grid.topology.coordinates;
    elements = grid.topology.cells;
    nElem = size(elements,1);
    stiffnesses = cell(nElem,1);

    for i = 1:nElem
        locCoords = coords(elements(i,:),:);
        bar = mean(locCoords, 1);
        dist = abs(dot(norm, bar) + d);

        blkID = int32(floor(dist / blockSize));

        D_eff = get_stiffness(blkID, E_rock, nu_rock, azimuths, dips, ...
            Kn_values, Ks_values, numbers_of_joints, blockSize);

        stiffnesses{i} = D_eff * 1e-6;   % save in MPa
    end
end
