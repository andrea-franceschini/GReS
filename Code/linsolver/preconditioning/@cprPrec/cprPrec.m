classdef cprPrec < preconditioner
% CPRPREC Constrained Pressure Residual (CPR) Preconditioner.
%
% This class implements a two-stage multiplicative Constrained Pressure
% Residual (CPR) preconditioner derived from the abstract 'preconditioner'
% base class.
%
% For coupled multiphase reservoir flow systems:
%   [ A_pp  A_ps ] [ x_p ]   [ b_p ]
%   [ A_sp  A_ss ] [ x_s ] = [ b_s ]
%
% Stage 1 (Elliptic / Pressure solve):
%   Decoupling to eliminate the primary saturation dependence from the
%   pressure equation:
%     - Quasi-IMPES ('quasi'):
%         alpha = - diag(A_ps) ./ diag(A_ss)
%     - True-IMPES ('true'):
%         alpha = rho_w(p) ./ rho_nw(p)  (volume-scaling from fluid PVT)
%     A_p = A_pp + diag(alpha) * A_sp
%     r_p = r_p + alpha .* r_s
%   The pressure correction x_p = A_p \ r_p is computed using Algebraic
%   Multigrid (aAMG).
%
% Stage 2 (Smoothing / Relaxation):
%   The intermediate residual is updated:
%     r^(1) = r - A * [x_p; 0]
%   and smoothed on the full system using ILU(0) (or other smoothers):
%     x^(2) = M_smooth \ r^(1)
%
% Multiplicative combination:
%     x = [x_p; 0] + x^(2)
%
% Parameter to switch Stage 1 decoupling:
%   obj.cprType = 'quasi' (Quasi-IMPES, default)
%   obj.cprType = 'true'  (True-IMPES)
%
% References:
%   Wallis, J. R. (1983). Incomplete Gaussian elimination as a preconditioning
%   for generalized conjugate gradient acceleration. SPE-12265-MS.
%   Cao, H. et al. (2005). Multi-stage preconditioners for parallel reservoir
%   simulation. SPE-93380-MS.

   properties (Access = public)

      % Decoupling strategy for Stage 1: 'quasi' (default) or 'true'
      cprType = 'quasi'

      % Smoother type: 'fsai' (default) or 'ilu'
      smootherType = 'fsai'

   end

   properties (Access = private)

      % Inner AMG for decoupled pressure equation
      AMGP = []

      % Inner FSAI smoother for full system
      FSAI = []

      % Factors for ILU(0) smoothing stage (alternative / fallback)
      L_ilu = []
      U_ilu = []

      % Problemsolver reference
      problemsolver

      % Discretizer domains
      domain
      
      % Number of domains/interfaces
      nDom
      nInt

      % Assembled system matrix used for residual updates
      Amat = []

   end

   properties (GetAccess = public, SetAccess = private)

      % Decoupling weights: alpha
      alpha = []

      % Decoupled pressure matrix: Ap = A11 + diag(alpha)*A21
      Ap = []

      % Number of pressure DOFs
      np = 0

      % Matrix symmetry (multiphase flow systems are nonsymmetric)
      PrecSym = false

      % AMG parameters structure
      params

      % Physics identifier (0.2 for multiphase flow)
      phys = 0.2

      % Max threads
      maxThreads

   end

   methods (Access = public)

      % Compute the preconditioner
      A = Compute(obj,A,varargin)

      % Left preconditioner application
      function x = ApplyLeft(obj,b,varargin)
         if isempty(obj.Amat)
            error('cprPrec must be computed before applying.');
         end

         np = obj.np;
         r1 = b(1:np);
         r2 = b(np+1:end);

         % 1. Restriction to decoupled pressure space: r_p = r_1 + alpha .* r_2
         rp = r1 + obj.alpha .* r2;

         % 2. Pressure solve with AMG (or direct solve fallback)
         if ~isempty(obj.AMGP) && ~isempty(obj.AMGP.Prec)
            xp = obj.AMGP.ApplyLeft(rp, obj.Ap);
         else
            xp = obj.Ap \ rp;
         end

         % 3. Prolongate pressure correction: x^(1) = [xp; 0]
         x1 = [xp; zeros(size(r2))];

         % 4. Residual update: r^(1) = b - A * x^(1)
         r_upd = b - obj.Amat * x1;

         % 5. Second-stage smoothing on updated residual
         switch lower(obj.smootherType)
            case {'fsai', 'afsai'}
               x2 = obj.FSAI.ApplyLeft(r_upd);
            case 'ilu'
               x2 = obj.U_ilu \ (obj.L_ilu \ r_upd);
            otherwise
               x2 = obj.FSAI.ApplyLeft(r_upd);
         end

         % 6. Multiplicative combination
         x = x1 + x2;
      end

      % Right preconditioner application (identity)
      function x = ApplyRight(~,b,varargin)
         x = b;
      end

      % Update matrix stored for residual update when full recomputation is skipped
      function updateStateBlocks(obj, A)
         if iscell(A)
            obj.Amat = cell2matrix(A);
         else
            obj.Amat = A;
         end
      end

      % Constructor
      function obj = cprPrec(debugflag,problemsolver,varargin)

         % Call base class constructor
         obj = obj@preconditioner();

         obj.DEBUGflag = debugflag;
         obj.problemsolver = problemsolver;
         obj.domain = problemsolver.domains;
         obj.nDom = problemsolver.nDom;
         obj.nInt = problemsolver.nInterf;

         % Check if user specified parameters in simparams
         if (isprop(problemsolver, 'simparams') || isfield(problemsolver, 'simparams')) && ...
            ~isempty(problemsolver.simparams) && ...
            (isprop(problemsolver.simparams, 'linSolverParams') || isfield(problemsolver.simparams, 'linSolverParams'))
            lsp = problemsolver.simparams.linSolverParams;
            if isstruct(lsp)
               if isfield(lsp, 'cprType')
                  obj.cprType = lsp.cprType;
               elseif isfield(lsp, 'decouplingType')
                  obj.cprType = lsp.decouplingType;
               end
               if isfield(lsp, 'smootherType')
                  obj.smootherType = lsp.smootherType;
               end
            end
         end

         % Parse varargin overrides
         if ~isempty(varargin)
            if isstruct(varargin{1})
               opts = varargin{1};
               if isfield(opts, 'cprType'), obj.cprType = opts.cprType; end
               if isfield(opts, 'smootherType'), obj.smootherType = opts.smootherType; end
            elseif ischar(varargin{1}) || isstring(varargin{1})
               for k = 1:2:numel(varargin)
                  param = lower(string(varargin{k}));
                  if k+1 <= numel(varargin)
                     val = varargin{k+1};
                     if ismember(param, ["cprtype", "type", "decouplingtype", "stage1", "stage1type"])
                        obj.cprType = val;
                     elseif ismember(param, ["smoothertype", "smoother", "stage2", "stage2type"])
                        obj.smootherType = val;
                     end
                  end
               end
            end
         end

         % Create inner aAMG instance for pressure
         obj.AMGP = aAMG(debugflag, problemsolver, "pressure");
         obj.maxThreads = obj.AMGP.maxThreads;
         obj.params = obj.AMGP.params;

         % Create inner aFSAI instance for smoothing
         obj.FSAI = aFSAI(debugflag);

         obj.PrecSym = false;

      end

   end

end
