function S_global = transform_joint_compliance(S_j, R)

    T = build_voigt_rotation_matrix(R);

    S_j_voigt = zeros(6);

    % Mapping from Marcio's code
    S_j_voigt(3,3) = S_j(3,3);  % normal compliance
    S_j_voigt(4,4) = S_j(2,2);  % shear
    S_j_voigt(5,5) = S_j(1,1);  % shear

    S_global = T' * S_j_voigt * T;
end
