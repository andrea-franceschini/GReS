function A = Compute(obj,A,varargin)
    % COMPUTE Assemble the two-stage Multiplicative CPR Preconditioner.
    %
    % Steps:
    % 1. Extract 2x2 blocks (pressure and saturation equations).
    % 2. Compute decoupling weights alpha based on obj.cprType:
    %    - 'quasi': Quasi-IMPES, alpha = -diag(A12) ./ diag(A22)
    %    - 'true':  True-IMPES,  alpha = rho_w(p) ./ rho_nw(p)
    % 3. Form decoupled pressure matrix: Ap = A11 + diag(alpha)*A21.
    % 4. Compute AMG hierarchy for the elliptic pressure matrix Ap.
    % 5. Compute FSAI of full system matrix Amat.
    % 6. Register Apply_L and Apply_R function handles.
    
    % Multiphase system is nonsymmetric
    obj.PrecSym = false;
    
    % Compute and apply Ruiz diagonal scaling if requested
    A = obj.Ruiz.Compute(A);
    
    % Extract block matrices and assembled sparse matrix
    A11 = A{1,1};
    A12 = A{1,2};
    A21 = A{2,1};
    A22 = A{2,2};
    obj.Amat = cell2matrix(A);
    
    ndofP = size(A11, 1);
    
    % Decoupling weights alpha and saturation factor d12 based on cprType
    % In GReS block ordering:
    % Row 1 = wetting mass balance, Row 2 = non-wetting mass balance
    % Col 1 = pressure (p), Col 2 = saturation (s)
    % Thus A12 = dr_w/ds, A22 = dr_nw/ds
    %      A11 = dr_w/dp, A21 = dr_nw/dp
    if strcmp(obj.cprType, 'true')
        % True-IMPES reduction via column-sum lumping (Eq 36a-36b, White et al. 2019):
        % Ass ≈ D_ss = diag(Ass^T * e), Aps ≈ D_ps = diag(Aps^T * e)
        e = ones(size(A12, 1), 1);
        d_ss = full(A12' * e);
        d_ps = full(A22' * e);
        alpha = zeros(size(d_ss));
        valid = abs(d_ss) > max(abs(d_ss))*eps;
        alpha(valid) = d_ps(valid) ./ d_ss(valid);
    
    else
        if ~strcmpi(obj.cprType,'quasi')
            warning('CPR:UnknownType', 'Unknown cprType "%s". Defaulting to "quasi".', obj.cprType);
        end

        % Quasi-IMPES reduction via diagonal extraction (Eq 37a-37b, White et al. 2019):
        % Ass ≈ D_ss = diagm(Ass), Aps ≈ D_ps = diagm(Aps)
        d_ss = full(diag(A12));
        d_ps = full(diag(A22));
        alpha = zeros(size(d_ss));
        valid = abs(d_ss) > max(abs(d_ss))*eps;
        alpha(valid) = d_ps(valid) ./ d_ss(valid);
    end
    
    obj.alpha = alpha(:);
    
    % Decoupled pressure operator (Eq 38, White et al. 2019):
    % S_pp = A21 - diag(alpha) * A11
    W = sparse(1:ndofP, 1:ndofP, obj.alpha, ndofP, ndofP);
    obj.Ap = A21 - W * A11;
    
    % Multiphase pressure operator with upwinded mobilities is non-symmetric
    symAp = false;
    
    % Test space: constant vector for pressure (scaled if Ruiz scaling is active)
    TV0 = ones(ndofP, 1);
    TV0 = obj.Ruiz.applyDinv(TV0, 1);
    
    % Stage 1: Compute AMG for decoupled pressure equation
    obj.AMGP.Compute(obj.Ap, symAp, TV0);
    
    % Stage 2: Compute smoother for full system
    obj.FSAI.Compute(obj.Amat, false);
    
    % Set application function handles
    obj.Apply_L = @(r) obj.ApplyLeft(r,ndofP,d_ss(:),A22);
    obj.Apply_R = @(r) obj.ApplyRight(r);

end
