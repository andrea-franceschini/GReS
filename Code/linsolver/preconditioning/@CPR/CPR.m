classdef CPR < preconditioner
% CPR Constrained Pressure Residual (CPR) Preconditioner.
%
% Derived from the abstract 'preconditioner' base class. Implements a
% two-stage multiplicative Constrained Pressure Residual (CPR) preconditioner
% for coupled multiphase flow systems.
%
% Stage 1 decouples the elliptic pressure equation by eliminating primary
% saturation dependence, then solves for the pressure update using an
% Algebraic Multigrid hierarchy (aAMG).
%
% Stage 2 updates the residual on the full coupled system and applies an
% approximate inverse smoother (aFSAI).
%
% The two corrections are combined multiplicatively to obtain the final
% preconditioned vector.

   properties (Access = public)

      % Decoupling strategy for Stage 1: 'quasi' (default) or 'true'
      cprType = 'quasi'

   end

   properties (SetAccess = public, GetAccess = public)

      % Inner AMG for decoupled pressure equation
      AMGP = []

      % Inner FSAI smoother for full system
      FSAI = []

   end

   properties (Access = private)

      % Problemsolver reference
      problemsolver

      % Discretizer domains
      domain

      % Number of domains/interfaces
      nDom
      nInt

   end

   properties (GetAccess = public, SetAccess = private)

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
      function x = ApplyLeft(obj,b,ndofP,d12,A22,Amat)

         r1 = b(1:ndofP);       % Wetting mass balance residual
         r2 = b(ndofP+1:end);   % Non-wetting mass balance residual

         % Jacobi saturation update
         zs = zeros(size(r1));
         valid = abs(d12) > eps*max(abs(d12));
         zs(valid) = r1(valid) ./ d12(valid);

         % Coupling update
         wp = r2 - A22 * zs;

         % AMG pressure solve
         zp = obj.AMGP.Apply_L(wp);

         % Intermediate solution and residual update
         x1 = [zp; zs];
         r_upd = b - Amat * x1;

         % Second-stage smoothing on updated residual using FSAI
         x2 = obj.FSAI.Apply_L(r_upd);

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

         % Check if user specified parameters in simparams
         if isfield(problemsolver.simparams, 'linSolverParams')
            lsp = problemsolver.simparams.linSolverParams;
            if isfield(lsp, 'cprType')
               obj.cprType = lsp.cprType;
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