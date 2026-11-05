function A = Compute(obj,A,varargin)
    % Multiphase system is nonsymmetric
    obj.PrecSym = false;
    
    % Compute and apply Ruiz diagonal scaling if requested
    A = obj.Ruiz.Compute(A);
    
    % Extract block matrices and assembled sparse matrix
    A11 = A{1,1};
    A12 = A{1,2};
    A21 = A{2,1};
    A22 = A{2,2};
    Amat = cell2matrix(A);
    ndofP = size(A11, 1);
    
    % Decoupling weights alpha and saturation factor d_ss based on cprType
    % Row 1 = wetting mass balance, Row 2 = non-wetting mass balance
    % Col 1 = pressure (p), Col 2 = saturation (s)
    if strcmp(obj.cprType, 'true')
        e = ones(size(A12, 1), 1);
        d_ss = full(A12' * e);
        d_ps = full(A22' * e);
        alpha = zeros(size(d_ss));
        valid = abs(d_ss) > max(abs(d_ss))*eps;
        alpha(valid) = d_ps(valid) ./ d_ss(valid);
    else
        d_ss = full(diag(A12));
        d_ps = full(diag(A22));
        alpha = zeros(size(d_ss));
        valid = abs(d_ss) > max(abs(d_ss))*eps;
        alpha(valid) = d_ps(valid) ./ d_ss(valid);
    end
    
    % Decoupled pressure operator:
    % S_pp = A21 - diag(alpha) * A11
    W = sparse(1:ndofP, 1:ndofP, alpha(:), ndofP, ndofP);
    Ap = A21 - W * A11;
    
    % Multiphase pressure operator with upwinded mobilities is non-symmetric
    symAp = false;
    
    % Test space: constant vector for pressure (scaled if Ruiz scaling is active)
    TV0 = ones(ndofP, 1);
    TV0 = obj.Ruiz.applyDinv(TV0, 1);
    
    % Compute AMG for decoupled pressure equation
    obj.AMGP.Compute(Ap, symAp, TV0);
    
    % Compute smoother for full system
    obj.FSAI.Compute(Amat, false);
    
    % Set application function handles
    obj.Apply_L = @(r) obj.ApplyLeft(r, ndofP, d_ss(:), A22, Amat);
    obj.Apply_R = @(r) obj.ApplyRight(r);

end
