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
% 5. Compute ILU(0) incomplete factorization of full system matrix Amat.
% 6. Register Apply_L and Apply_R function handles.

   % Multiphase system is nonsymmetric
   obj.PrecSym = false;

   % Compute and apply Ruiz diagonal scaling if requested
   if ~isempty(obj.Ruiz) && isa(obj.Ruiz, 'RuizScaling') && obj.Ruiz.scalingFlag
      A = obj.Ruiz.Compute(A);
   end

   % Extract block matrices and assembled sparse matrix
   if iscell(A)
      A11 = A{1,1};
      A12 = A{1,2};
      A21 = A{2,1};
      A22 = A{2,2};
      obj.Amat = cell2matrix(A);
   else
      obj.Amat = A;
      np = obj.domain(1).dofm.getNumbDoF("pressure");
      A11 = A(1:np, 1:np);
      A12 = A(1:np, np+1:end);
      A21 = A(np+1:end, 1:np);
      A22 = A(np+1:end, np+1:end);
   end

   obj.np = size(A11, 1);
   np = obj.np;

   % Decoupling weights alpha and saturation factor d12 based on cprType
   % In GReS block ordering:
   % Row 1 = wetting mass balance, Row 2 = non-wetting mass balance
   % Col 1 = pressure (p), Col 2 = saturation (s)
   % Thus A12 = dr_w/ds (Ass in White et al.), A22 = dr_nw/ds (Aps in White et al.)
   %      A11 = dr_w/dp (Asp in White et al.), A21 = dr_nw/dp (App in White et al.)
   switch lower(string(obj.cprType))
      case {"true", "true-impes", "trueimpes"}
         % True-IMPES reduction via column-sum lumping (Eq 36a-36b, White et al. 2019):
         % Ass ≈ D_ss = diag(Ass^T * e), Aps ≈ D_ps = diag(Aps^T * e)
         e = ones(size(A12, 1), 1);
         d_ss = full(A12' * e);
         d_ps = full(A22' * e);
         alpha = zeros(size(d_ss));
         valid = abs(d_ss) > 1e-14;
         alpha(valid) = d_ps(valid) ./ d_ss(valid);

      case {"quasi", "quasi-impes", "quasiimpes"}
         % Quasi-IMPES reduction via diagonal extraction (Eq 37a-37b, White et al. 2019):
         % Ass ≈ D_ss = diagm(Ass), Aps ≈ D_ps = diagm(Aps)
         d_ss = full(diag(A12));
         d_ps = full(diag(A22));
         alpha = zeros(size(d_ss));
         valid = abs(d_ss) > 1e-14;
         alpha(valid) = d_ps(valid) ./ d_ss(valid);

      case {"pvt", "physical"}
         % Physical volume-scaling using fluid phase densities rho_nw(p) / rho_w(p)
         alpha = [];
         if ~isempty(obj.domain)
            try
               tpSolver = obj.domain(1).getPhysicsSolver('TwoPhaseFlow');
               if ~isempty(tpSolver) && ~isempty(tpSolver.props.wet) && ~isempty(tpSolver.props.nwet)
                  p = obj.domain(1).getState("pressure");
                  cIdx = obj.domain(1).dofm.getFieldCells("pressure");
                  pCells = p(cIdx);
                  rhoW = tpSolver.props.wet.density(pCells);
                  rhoNW = tpSolver.props.nwet.density(pCells);
                  alpha = - full(value(rhoNW ./ rhoW));
               end
            catch
               alpha = [];
            end
         end
         d_ss = full(diag(A12));
         if isempty(alpha)
            d_ps = full(diag(A22));
            alpha = zeros(size(d_ss));
            valid = abs(d_ss) > 1e-14;
            alpha(valid) = d_ps(valid) ./ d_ss(valid);
         end

      otherwise
         warning('CPR:UnknownType', 'Unknown cprType "%s". Defaulting to "quasi".', obj.cprType);
         d_ss = full(diag(A12));
         d_ps = full(diag(A22));
         alpha = zeros(size(d_ss));
         valid = abs(d_ss) > 1e-14;
         alpha(valid) = d_ps(valid) ./ d_ss(valid);
   end

   obj.alpha = alpha(:);
   obj.d12 = d_ss(:);
   obj.A22 = A22;

   % Decoupled pressure operator (Eq 38, White et al. 2019):
   % S_pp = A21 - diag(alpha) * A11
   W = sparse(1:np, 1:np, obj.alpha, np, np);
   obj.Ap = A21 - W * A11;

   % Multiphase pressure operator with upwinded mobilities is non-symmetric
   symAp = false;

   % Test space: constant vector for pressure (scaled if Ruiz scaling is active)
   TV0 = ones(np, 1);
   if ~isempty(obj.Ruiz) && isa(obj.Ruiz, 'RuizScaling') && obj.Ruiz.scalingFlag
      TV0 = obj.Ruiz.applyDinv(TV0, 1);
   end

   % Stage 1: Compute AMG for decoupled pressure equation
   if ~isempty(obj.AMGP)
      obj.AMGP.Compute(obj.Ap, symAp, TV0);
   end

   % Stage 2: Compute smoother for full system (FSAI or ILU)
   switch lower(obj.smootherType)
      case {'fsai', 'afsai'}
         if isempty(obj.FSAI)
            obj.FSAI = aFSAI(obj.DEBUGflag);
         end
         obj.FSAI.Compute(obj.Amat, false);

      case 'ilu'
         setup = struct('type', 'nofill');
         try
            [obj.L_ilu, obj.U_ilu] = ilu(obj.Amat, setup);
         catch
            setup.udiag = 1;
            [obj.L_ilu, obj.U_ilu] = ilu(obj.Amat, setup);
         end

      otherwise
         if isempty(obj.FSAI)
            obj.FSAI = aFSAI(obj.DEBUGflag);
         end
         obj.FSAI.Compute(obj.Amat, false);
   end

   % Set application function handles
   obj.Apply_L = @(r) obj.ApplyLeft(r);
   obj.Apply_R = @(r) obj.ApplyRight(r);

end
