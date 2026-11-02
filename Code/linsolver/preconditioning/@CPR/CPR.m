classdef CPR < preconditioner
% CPR Constrained Pressure Residual (CPR) Preconditioner.
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
%   and smoothed on the full system using an adaptive FSAI:
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

   end

   properties (SetAccess = public, GetAccess = public)

      % Inner AMG for decoupled pressure equation
      AMGP = []

      % Inner FSAI smoother for full system
      FSAI = []

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

      % Matrix symmetry (multiphase flow systems are nonsymmetric)
      PrecSym = false

      % AMG parameters structure
      params

      % Physics
      phys

      % Max threads
      maxThreads

   end

   methods (Access = public)

      % Compute the preconditioner
      A = Compute(obj,A,varargin)

      % Left preconditioner application
      function x = ApplyLeft(obj,b,varargin)

         % Number of degrees of freedom of the pressure
         ndofP = varargin{1};

         % Saturation diagonal / column-sum factor
         d12 = varargin{2};

         % Coupling block: A22
         A22 = varargin{3};

         r1 = b(1:ndofP);       % Wetting mass balance residual
         r2 = b(ndofP+1:end);   % Non-wetting mass balance residual

         % Jacobi saturation update zs = (D_ss^{(CPR)})^{-1} * r1
         zs = zeros(size(r1));
         valid_d12 = abs(d12) > eps*max(abs(d12));
         zs(valid_d12) = r1(valid_d12) ./ d12(valid_d12);

         % Coupling update wp = r2 - A22 * zs
         wp = r2 - A22 * zs;

         % AMG pressure solve: Ap * zp = wp
         zp = obj.AMGP.ApplyLeft(wp, obj.Ap);

         % Prolongate both pressure and saturation: x^(1) = [zp; zs]
         x1 = [zp; zs];

         % Residual update r^(1) = b - Amat * x^(1)
         r_upd = b - obj.Amat * x1;

         % Second-stage smoothing on updated residual using FSAI
         x2 = obj.FSAI.ApplyLeft(r_upd);

         % Multiplicative combination
         x = x1 + x2;
      end

      % Right preconditioner application (identity)
      function x = ApplyRight(~,b,varargin)
         x = b;
      end

      % Constructor
      function obj = CPR(debugflag,problemsolver,varargin)

         % Call base class constructor
         obj = obj@preconditioner();

         obj.DEBUGflag = debugflag;
         obj.problemsolver = problemsolver;
         obj.domain = problemsolver.domains;
         obj.nDom = problemsolver.nDom;
         obj.nInt = problemsolver.nInterf;

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
