function D_eff = get_stiffness(blkID, E_rock, nu_rock, ...
    azimuths, dips, Kn, Ks, numbers_of_joints, block_size)

    % Intact rock compliance
    S_r = compliance_intact_rock(E_rock, nu_rock);

    %tot_num_joints = 0;
    %max_num_joints = 0;

    if blkID < size(numbers_of_joints,1)

        S_j_total = zeros(6);

        for i = 1:length(azimuths.mean)

            % Randomize orientation
            azimuth = generate_truncated_random_numbers( ...
                azimuths.mean(i), azimuths.std(i), 0, inf);

            dip = generate_truncated_random_numbers( ...
                dips.mean(i), dips.std(i), 0, inf);

            % Truncated stiffnesses
            Kn_value = generate_truncated_random_numbers( ...
                Kn.mean(i), Kn.std(i), Kn.min(i), Kn.max(i));

            Ks_value = generate_truncated_random_numbers( ...
                Ks.mean(i), Ks.std(i), Ks.min(i), Ks.max(i));

            % Number of joints
            mu_joints = numbers_of_joints(blkID+1, i);
            num_joints = round(generate_truncated_random_numbers( ...
                mu_joints, 1, 0, inf));

            if num_joints > 0

                spacing = block_size / num_joints;

                [alpha, beta, gamma] = ...
                    az_dip_to_euler_xyz(azimuth, dip);

                R = rotation_matrix(alpha, beta, gamma);

                S_j_local = joint_compliance(Kn_value, Ks_value);

                S_j_global = ...
                    transform_joint_compliance(S_j_local, R) ...
                    / spacing;

                S_j_total = S_j_total + S_j_global;
            end

            %tot_num_joints = tot_num_joints + num_joints;
            %max_num_joints = max(max_num_joints, num_joints);

        end

        S_eff = S_r + S_j_total;
        D_eff = inv(S_eff);

    else
        D_eff = inv(S_r);
    end
end
