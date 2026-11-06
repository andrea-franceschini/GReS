classdef fixedStress < preconditioner

   properties (SetAccess = protected, GetAccess = public)

      % Inner preconditioners
      AMGMech = []
      FlowPrec = []

      % Flow preconditioner type: 'AMG' (single-phase) or 'CPR' (two-phase)
      flowPrecType = 'AMG'

      % Flag indicating if mechanics is linear (allows AMGMech recycling)
      isMechLinear = false

      % Flow part setup time (excluding one-time linear mechanics AMG)
      tSetupFlow = 0

   end

   properties (Access = private)

      % Problemsolver params
      problemsolver

      % Discretizer
      domain
      
      % Number of domains/interfaces
      nDom
      nInt

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

      % Getter for the function handle to apply the left preconditioner
      function x = ApplyLeft(obj,b,varargin)
         if strcmpi(obj.flowPrecType,'CPR')
            n1 = varargin{1};
            Afu = varargin{2};

            bu = b(1:n1);
            bf = b(n1+1:end);

            % White et al. (2019) block triangular solve:
            % Mechanics solve: zu = M_uu^{-1} bu
            zu = obj.AMGMech.Apply_L(bu);

            % Flow coupling update: yf = bf - Afu * zu
            yf = bf - Afu * zu;

            % Flow CPR solve: zf = CPR(yf)
            zf = obj.FlowPrec.Apply_L(yf);

            % Combine displacement and flow solutions
            x = [zu; zf];

         else
            A11 = varargin{1};
            B1 = varargin{2};

            % Get mechanics size
            n1 = size(A11,1);

            % Apply AMG for fluid part (approximate Schur complement)
            x2 = obj.FlowPrec.Apply_L(b(n1+1:end));

            % Apply the top part of the block preconditioner (upper block triangular solve)
            x11 = b(1:n1) - B1*x2;
            x1 = obj.AMGMech.Apply_L(x11);

            % Compose the solution
            x = [x1;x2];
         end
      end

      % Getter for the function handle to apply the right preconditioner
      function x = ApplyRight(~,b,varargin)
         if nargin < 2
            error('Not enough arguments for fixedStress apply Right');
         end

         x = b;
      end

      % Return setup time for recomputing the preconditioner
      % (for linear mechanics, only the flow part is recomputed, not the mechanics part)
      function t = getRecomputeSetupTime(obj, defaultTSetup)
         if obj.isMechLinear && ~isempty(obj.tSetupFlow) && obj.tSetupFlow > 0
            t = obj.tSetupFlow;
         else
            t = defaultTSetup;
         end
      end

      % Constructor Function
      function obj = fixedStress(debugflag,problemsolver,varargin)

         % Call the constructor of the abstract class
         obj = obj@preconditioner();
         
         % Set the debugflag
         obj.DEBUGflag = debugflag;

         % Get the domains
         obj.problemsolver = problemsolver;
         obj.domain = problemsolver.domains;

         obj.nInt = problemsolver.nInterf;
         obj.nDom = problemsolver.nDom;

         % Check if two-phase flow is present
         hasSat = false;
         for d = 1:numel(obj.domain)
            if ismember("saturation", obj.domain(d).dofm.getVariableNames())
               hasSat = true;
               break;
            end
         end

         % Check if mechanics is linear (allows AMGMech recycling)
         obj.isMechLinear = obj.checkMechLinear();

         % Create the inner mechanics AMG preconditioner
         obj.AMGMech = aAMG(debugflag,problemsolver,"displacements");

         if hasSat
            obj.flowPrecType = 'CPR';
            obj.FlowPrec = CPR(debugflag, problemsolver);
            obj.maxThreads = obj.FlowPrec.maxThreads;
            obj.params = obj.FlowPrec.params;
         else
            obj.flowPrecType = 'AMG';
            obj.FlowPrec = aAMG(debugflag,problemsolver,"pressure");
            obj.maxThreads = obj.FlowPrec.maxThreads;
            obj.params = obj.FlowPrec.params;
         end
      end
   end

   methods (Access = private)
      function isLin = checkMechLinear(obj)
         isLin = false;
         for d = 1:numel(obj.domain)
            dom = obj.domain(d);
            if ismember("MultiphasePoromechanics", dom.solverNames)
               mSolver = dom.getPhysicsSolver("MultiphasePoromechanics");
               isLin = mSolver.getMechSolver().isLinear();
               return;
            elseif ismember("BiotFullyCoupled", dom.solverNames)
               mSolver = dom.getPhysicsSolver("BiotFullyCoupled");
               isLin = mSolver.getMechSolver().isLinear();
               return;
            end
         end
      end
   end
end
