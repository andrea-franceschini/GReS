classdef fixedStress < preconditioner

   properties (SetAccess = protected, GetAccess = public)

      % Inner preconditioners
      AMGMech = []
      AMGFlux = []
      FlowPrec = []

      % Flow preconditioner type: 'AMG' (single-phase) or 'CPR' (two-phase)
      flowPrecType = 'AMG'

      % Cached blocks for coupled solve
      Auu = []
      Auf = []
      Afu = []
      S_ff = []

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
         switch upper(obj.flowPrecType)
            case 'CPR'
               n1 = size(obj.Auu, 1);
               bu = b(1:n1);
               bf = b(n1+1:end);

               % Algorithm (White et al. 2019):
               % Mechanics solve: zu = M_uu^{-1} bu
               zu = obj.AMGMech.ApplyLeft(bu, obj.Auu);

               % Flow coupling update: yf = bf - Afu * zu
               yf = bf - obj.Afu * zu;

               % Flow CPR solve: zf = CPR(yf)
               zf = obj.FlowPrec.ApplyLeft(yf);

               % Combine
               x = [zu; zf];

            otherwise
               if nargin >= 4
                  S = varargin{1};
                  A11 = varargin{2};
                  B1 = varargin{3};
               elseif ~isempty(obj.S_ff) && ~isempty(obj.Auu) && ~isempty(obj.Auf)
                  S = obj.S_ff;
                  A11 = obj.Auu;
                  B1 = obj.Auf;
               else
                  error('Not enough arguments for fixedStress apply Left');
               end

               % Get mechanics size
               n1 = size(A11,1);

               % Apply AMG for fluid part
               x2 = obj.AMGFlux.ApplyLeft(b(n1+1:end),S);

               % Apply the top part of the block preconditioner
               x11 = b(1:n1) - B1*x2;
               x1 = obj.AMGMech.ApplyLeft(x11,A11);

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

      % % Update state blocks when preconditioner recomputation is skipped
      % function updateStateBlocks(obj, A)
      %    if strcmpi(obj.flowPrecType, 'CPR')
      %       if iscell(A)
      %          obj.Auu = A{1,1};
      %          obj.Auf = [A{1,2}, A{1,3}];
      %          obj.Afu = [A{2,1}; A{3,1}];
      %          Aff = {A{2,2}, A{2,3}; A{3,2}, A{3,3}};
      %          if ~isempty(obj.FlowPrec)
      %             obj.FlowPrec.updateStateBlocks(Aff);
      %          end
      %       end
      %    else
      %       if iscell(A)
      %          obj.Auu = A{1,1};
      %          obj.Auf = A{1,2};
      %       end
      %    end
      % end

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
            obj.AMGFlux = aAMG(debugflag,problemsolver,"pressure");
            obj.maxThreads = obj.AMGFlux.maxThreads;
            obj.params = obj.AMGFlux.params;
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
            elseif ismember("Poromechanics", dom.solverNames)
               mSolver = dom.getPhysicsSolver("Poromechanics");
               isLin = mSolver.isLinear();
               return;
            elseif ismember("EmbeddedFractureMechanics", dom.solverNames)
               mSolver = dom.getPhysicsSolver("EmbeddedFractureMechanics");
               isLin = mSolver.isLinear();
               return;
            end
         end
      end
   end
end
