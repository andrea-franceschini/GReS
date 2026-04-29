function D_eff = get_stiffness(blkID, E_rock, nu_rock, ...
    azimuths, dips, Kn, Ks, numbers_of_joints, block_size)

    % Intact rock compliance
    S_r = compliance_intact_rock(E_rock, nu_rock);

    if blkID < size(numbers_of_joints,1)

        S_j_total = zeros(6);

        for i = 1:length(azimuths)

            azimuth = azimuths(i);
            dip = dips(i);
            Kn_value = Kn(i);
            Ks_value = Ks(i);

            % Number of joints
            num_joints = numbers_of_joints(blkID+1, i);

            if num_joints > 0

                spacing = block_size / num_joints;

                [alpha, beta, gamma] = ...
                    az_dip_to_euler_xyz(azimuth, dip);

                R = rotation_matrix(alpha, beta, gamma);

                S_j_local = joint_compliance(Kn_value, Ks_value);

                S_j_global = ...
                    transform_joint_compliance(S_j_local, R) / spacing;

                S_j_total = S_j_total + S_j_global;
            end

        end

        S_eff = S_r + S_j_total;
        D_eff = inv(S_eff);

    else
        D_eff = inv(S_r);
    end
end
