function gradN = computeShapeFunctionGradients(coords)
% Computes ∇N for a linear tetrahedral element.
% Input:  coords - 4x3 matrix [x y z] for the 4 nodes
% Output: gradN  - 4x3 matrix, each row = [dNi/dx dNi/dy dNi/dz]

    x = coords(:,1);
    y = coords(:,2);
    z = coords(:,3);

    % Volume (6 * V = det(M))
    M = [ones(4,1), x, y, z];
    V6 = det(M);  % = 6*V
    if abs(V6) < 1e-12
        error('Degenerate element with near-zero volume.');
    end

    gradN = zeros(4,3);
    for i = 1:4
        % Delete i-th row of M
        Mi = M;
        Mi(i,:) = [];
        % Cofactor terms for x,y,z
        dN_dx =  det(Mi(:,[2 3 4]));
        dN_dy = -det(Mi(:,[1 3 4]));
        dN_dz =  det(Mi(:,[1 2 4]));
        gradN(i,:) = (-1)^(i+1) * [dN_dx, dN_dy, dN_dz] / V6;
    end
end