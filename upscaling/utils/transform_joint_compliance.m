function S_global = transform_joint_compliance(S_j, R)

    T = voigt_rotation_matrix(R);

    S_j_voigt = zeros(6);

    % Mapping from Marcio's code
    % the joint's normal is aligned with x, thus:
    % n -> x
    % t -> xz
    % t -> xy
    % it will be rotated by T
    S_j_voigt(1,1) = S_j(3,3);  % normal compliance
    S_j_voigt(5,5) = S_j(2,2);  % shear
    S_j_voigt(6,6) = S_j(1,1);  % shear

    S_global = T' * S_j_voigt * T;
end
