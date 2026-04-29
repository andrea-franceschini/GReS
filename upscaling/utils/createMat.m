function stiffnesses = createMat(grid, meshProp, rock, blockSize, jointFamilies)

    % Get joint/rock parameters
    [E_rock, nu_rock, azimuths, dips, Kn_values, Ks_values, ...
        numbers_of_joints] = set_parameters(rock, jointFamilies);

    theta = deg2rad(meshProp.fault_angle);
    cos_t = cos(theta);
    sin_t = sin(theta);
    pt   = meshProp.center;
    norm = [cos_t, 0, sin_t];
    d    = -dot(norm, pt);

    coords = grid.coordinates;
    elements = grid.cells.connectivity;
    nElem = size(elements,1);

    azimuthVec = zeros(nElem,length(azimuths.mean));
    dipVec = zeros(nElem,length(azimuths.mean));
    Kn_valueVec = zeros(nElem,length(azimuths.mean));
    Ks_valueVec = zeros(nElem,length(azimuths.mean));
    for i = 1:length(azimuths.mean)
        % Randomize orientation
        azimuthVec(:,i) = truncated_random( ...
            azimuths.mean(i), azimuths.std(i), 0, inf, nElem);

        dipVec(:,i) = truncated_random( ...
            dips.mean(i), dips.std(i), 0, inf, nElem);

        % Truncated stiffnesses
        Kn_valueVec(:,i) = truncated_random( ...
            Kn_values.mean(i), Kn_values.std(i), Kn_values.min(i), ...
            Kn_values.max(i), nElem);

        Ks_valueVec(:,i) = truncated_random( ...
            Ks_values.mean(i), Ks_values.std(i), Ks_values.min(i), ...
            Ks_values.max(i), nElem);
    end

    stiffnesses = cell(nElem,1);

    for i = 1:nElem
        locCoords = coords(elements(i,:),:);
        bar = mean(locCoords, 1);
        dist = abs(dot(norm, bar) + d);

        blkID = floor(dist / blockSize);

        D_eff = get_stiffness(blkID, E_rock, nu_rock, azimuthVec(i,:), ...
            dipVec(i,:), Kn_valueVec(i,:), Ks_valueVec(i,:), ...
            numbers_of_joints, blockSize);

        stiffnesses{i} = D_eff * 1e-6;   % save in MPa
    end
end
