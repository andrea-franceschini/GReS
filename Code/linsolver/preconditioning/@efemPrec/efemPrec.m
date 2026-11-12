classdef efemPrec < preconditioner

   properties (SetAccess = private,GetAccess=public)

      % Inner preconditioner
      AMG = []

      % A22 block Inverter
      FSAI = []

      % Problemsolver params
      problemsolver

      % Discretizer
      domain
      
      % Number of domains/interfaces
      nDom
      nInt

      % Schur complement matrix on which AMG was computed on
      S

      % Nonsymmetric matrix tolerance
      nsyTol

      % Linear solver tolerance to compute the correct regularization
      linsolverTol

      % FSAI params
      nstep = 15
   end

   properties (Access = public)
      % Regularization method for A22: 'diagShift' (default) or 'racp'
      % Note: Change flag directly here on file (no XML user change allowed)
      regularizationMethod = 'diagShift'

      % RACP augmentation parameter (default: 1e-3)
      omegaRACP = 1e-3
   end

   properties (GetAccess = public,SetAccess = private)
      % Symmetry of the matrix on which the preconditioner has been
      % computed
      PrecSym = true

      % AMG params structure
      params

      % Physics
      phys

      % Max Threads
      maxThreads
   end

   methods
      % Function to compute the preconditioner
      A = Compute(obj,A,sym,varargin)

      % Set diagonal block scope and propagate to AMG and FSAI children
      function setBlockScope(obj, sysInfo, blockIds)
         setBlockScope@preconditioner(obj, sysInfo, blockIds);
         uBlk = obj.getVarBlocks("displacements");
         jBlk = obj.getVarBlocks("fractureJump");
         obj.AMG.setBlockScope(sysInfo, obj.precBlocks(uBlk));
         obj.FSAI.setBlockScope(sysInfo, obj.precBlocks(jBlk));
      end

      % Getter for the function handle to apply the left preconditioner
      function x = ApplyLeft(obj,b,varargin)
         if nargin < 4
            error('Not enough arguments for fixedStress apply Left');
         end

         B1 = varargin{1};
         B2 = varargin{2};
         invC = varargin{3};
         
         % Get mechanics size
         n1 = size(obj.S,1);

         % Apply FSAI for the 22 block solution
         b2 = b(n1+1:end);
         x2 = invC(b2);

         % Apply the coupling
         x11 = b(1:n1) - B1*x2;
         
         % Apply AMG to complete the block upper Gauss Seidel
         x1 = obj.AMG.ApplyLeft(x11,obj.S);

         % Apply the coupling correction to the state block
         x2 = invC(b2 - B2*x1);

         % Compose the solution
         x = [x1;x2];
         
      end

      % Getter for the function handle to apply the right preconditioner
      function x = ApplyRight(obj,b,varargin)
         if nargin < 2
            error('Not enough arguments for fixedStress apply Right');
         end

         x = b;
      end

      % Make A22 slightly more diagonally dominant
      function A22 = regularizeA22(obj,A22)
         d = abs(diag(A22));
         maxD = max(d);
         if maxD == 0, maxD = 1; end
         tol = obj.linsolverTol;
         
         diagShift = (2.2e-15 / (tol * maxD)) * d;
         A22 = A22 + diag(diagShift);
      end

      % Compute physical stiffness scaling k_raw on unscaled blocks
      function k_raw = computePhysicalK(~, A11, A12)
         nw = size(A12, 2);
         d11 = full(diag(A11));
         mean_d11 = mean(d11);
         if mean_d11 == 0, mean_d11 = 1.0; end
         k_raw = zeros(nw, 1);
         for j = 1:nw
            v = A12(:, j);
            supp = find(v);
            if ~isempty(supp)
               m_a = max(d11(supp));
               if m_a <= 0, m_a = mean_d11; end
               m_b = full(v(supp)' * v(supp));
               k_raw(j) = m_b / m_a;
            else
               k_raw(j) = 1.0;
            end
         end
      end

      % Physics-based RACP augmentation for A22 block computed directly on Ruiz-scaled blocks
      function A22 = augmentA22_RACP(obj, A11, A12, A21, A22, k_phys)
         if nargin < 6 || isempty(k_phys)
            % If not provided, compute from Ruiz blocks or unscaled
            d11 = full(diag(A11));
            mean_d11 = mean(d11);
            if mean_d11 == 0, mean_d11 = 1.0; end
            nw = size(A22, 1);
            k_phys = zeros(nw, 1);
            for j = 1:nw
               v = A12(:, j);
               supp = find(v);
               if ~isempty(supp)
                  m_a = max(d11(supp));
                  if m_a <= 0, m_a = mean_d11; end
                  m_b = full(v(supp)' * v(supp));
                  k_phys(j) = m_b / m_a;
               else
                  k_phys(j) = 1.0;
               end
            end
         end

         % Compute relative coupling normalized to [0, 1]
         maxK = max(k_phys);
         if maxK <= 0, maxK = 1.0; end
         k_rel = k_phys / maxK;

         % Base diagShift magnitude in Ruiz space
         d = abs(diag(A22));
         maxD = max(d);
         if maxD == 0, maxD = 1; end
         tol = obj.linsolverTol;

         % RACP-weighted shift: preserves exact stable numerical scale of diagShift
         % while scaling the shift according to local rock-fracture coupling energy
         diagShift = (2.2e-15 / (tol * maxD)) * (k_rel .* d);
         A22 = A22 + diag(sparse(diagShift));
      end

      % Dispatches between diagShift and racp based on regularizationMethod
      function A22 = regularizeOrAugmentA22(obj, A11, A12, A21, A22)
         switch lower(obj.regularizationMethod)
            case 'racp'
               A22 = obj.augmentA22_RACP(A11, A12, A21, A22);
            case 'diagshift'
               A22 = obj.regularizeA22(A22);
            otherwise
               error('Unknown regularizationMethod: %s. Use ''diagShift'' or ''racp''.', obj.regularizationMethod);
         end
      end
      
      function updateStateBlocks(obj, A)
         % Cheap update for reuse

         % Check symmetry
         A22 = A{2,2};
         symm = norm(A22-A22','f')/norm(A22,'f') < obj.nsyTol;

         % Make A22 slightly more diagonally dominant or augment with RACP
         A22 = obj.regularizeOrAugmentA22(A{1,1}, A{1,2}, A{2,1}, A{2,2});

         % Compute the factorization of the 22 block
         obj.FSAI.Compute(A22,symm,true);
      
         % Compute reverse Schur complement of A11
         invC = @(x) obj.FSAI.ApplyLeft(x);
         B1 = A{1,2};
         B2 = A{2,1};
      
         % Update function handle while keeping the frozen AMG hierarchy S
         obj.Apply_L = @(x) obj.ApplyLeft(x, B1, B2, invC);
      end

      % Constructor Function
      function obj = efemPrec(debugflag,problemsolver,nsyTol)

         % Call the constructor of the abstract class
         obj = obj@preconditioner();
         
         % Set the debugflag
         obj.DEBUGflag = debugflag;

         % Get the domains
         obj.problemsolver = problemsolver;
         obj.domain = problemsolver.domains;

         obj.nInt = problemsolver.nInterf;
         obj.nDom = problemsolver.nDom;

         % Create the inner AMG preconditioner
         obj.AMG = aAMG(debugflag,problemsolver,"displacements");

         % Create the inner FSAI solver
         obj.FSAI = aFSAI(debugflag);

         % Get the parameters
         obj.maxThreads = obj.AMG.maxThreads;
         obj.params = obj.AMG.params;
         obj.nsyTol = nsyTol;
         obj.linsolverTol = problemsolver.simparams.relTol;
      end
   end
end
