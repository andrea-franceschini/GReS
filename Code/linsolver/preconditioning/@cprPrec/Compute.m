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

   % Decoupling weights alpha based on cprType
   switch lower(string(obj.cprType))
      case {"true", "true-impes", "trueimpes"}
         % True-IMPES: physical volume-scaling using fluid phase densities
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
                  alpha = full(value(rhoW ./ rhoNW));
               end
            catch
               alpha = [];
            end
         end

         % Fallback to Quasi-IMPES if fluid properties unavailable
         if isempty(alpha)
            d12 = full(diag(A12));
            d22 = full(diag(A22));
            alpha = zeros(size(d22));
            valid = abs(d22) > 1e-14;
            alpha(valid) = -d12(valid) ./ d22(valid);
         end

      case {"quasi", "quasi-impes", "quasiimpes"}
         % Quasi-IMPES: algebraic diagonal decoupling from Jacobian
         d12 = full(diag(A12));
         d22 = full(diag(A22));
         alpha = zeros(size(d22));
         valid = abs(d22) > 1e-14;
         alpha(valid) = -d12(valid) ./ d22(valid);

      otherwise
         warning('cprPrec:UnknownType', 'Unknown cprType "%s". Defaulting to "quasi".', obj.cprType);
         d12 = full(diag(A12));
         d22 = full(diag(A22));
         alpha = zeros(size(d22));
         valid = abs(d22) > 1e-14;
         alpha(valid) = -d12(valid) ./ d22(valid);
   end

   obj.alpha = alpha(:);

   % Decoupled pressure operator: Ap = A11 + diag(alpha)*A21
   W = sparse(1:np, 1:np, obj.alpha, np, np);
   obj.Ap = A11 + W * A21;

   % Determine symmetry of Ap
   normAp = norm(obj.Ap, 'fro');
   if normAp > 1e-12
      symAp = (norm(obj.Ap - obj.Ap', 'fro') / normAp) < 0.05;
   else
      symAp = false;
   end

   % Test space: constant vector for pressure
   TV0 = ones(np, 1);

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
