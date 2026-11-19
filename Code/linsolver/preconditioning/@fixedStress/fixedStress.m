classdef fixedStress < preconditioner

   properties (SetAccess = protected, GetAccess = public)

      % Inner preconditioners
      MechPrec = []
      FlowPrec = []

      % Flow preconditioner type: 'AMG' (single-phase) or 'CPR' (two-phase)
      flowPrecType = 'AMG'

      % Flow part setup time (excluding one-time linear mechanics)
      tSetupFlow = 0

      % Check if the old preconditioner was computed for linear physics
      isMechLinearOld = false;

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
         error('Undefined, use the specific functions ApplyLeft(...)');
      end

      function x = ApplyLeftCPR(obj,b,varargin)
         n1 = varargin{1};
         Afu = varargin{2};

         bu = b(1:n1);
         bf = b(n1+1:end);

         % White et al. (2019) block triangular solve:
         % Mechanics solve: zu = M_uu^{-1} bu
         zu = obj.MechPrec.Apply_L(bu);

         % Flow coupling update: yf = bf - Afu * zu
         yf = bf - Afu * zu;

         % Flow CPR solve: zf = CPR(yf)
         zf = obj.FlowPrec.Apply_L(yf);

         % Combine displacement and flow solutions
         x = [zu; zf];
      end


      function x = ApplyLeftAMG(obj,b,varargin)
         A11 = varargin{1};
         B1 = varargin{2};

         % Get mechanics size
         n1 = size(A11,1);

         % Apply AMG for fluid part (approximate Schur complement)
         x2 = obj.FlowPrec.Apply_L(b(n1+1:end));

         % Apply the top part of the block preconditioner (upper block triangular solve)
         x11 = b(1:n1) - B1*x2;
         x1 = obj.MechPrec.Apply_L(x11);

         % Compose the solution
         x = [x1;x2];
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

         % Create the inner mechanics AMG preconditioner
         obj.MechPrec = aAMG(debugflag,problemsolver,"displacements");

         % The two versions supported for single domain are the following,
         % but RACP as flowPrec is supported as well in the case of
         % multiple domains
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

      % Set sub-preconditioners for mechanics and flow
      function setSubPrecs(obj, varargin)
         for k = 1:2:numel(varargin)
            key = string(varargin{k});
            if strcmpi(key, "mechanics")
               obj.MechPrec = varargin{k+1};
            elseif strcmpi(key, "flow")
               obj.FlowPrec = varargin{k+1};
               if isa(obj.FlowPrec, 'CPR')
                  obj.flowPrecType = 'CPR';
               elseif isa(obj.FlowPrec, 'aAMG')
                  obj.flowPrecType = 'AMG';
               elseif isa(obj.FlowPrec, 'RACP')
                  obj.flowPrecType = 'RACP';
               else
                  error('Unsupported flow preconditioner');
               end
            end
         end
      end

      % Determine adaptive block ordering for activeBlockIds
      function blockOrder = getBlockOrdering(obj, blocks, activeBlockIds)
         activeMeta = blocks(activeBlockIds);
         
         % Get the flow blocks (the interfaces in varName have the
         % variables that they are coupling)
         isFlow = ismember([activeMeta.varName], ["pressure", "saturation"]);
         
         % Give the flow blocks to the flow preconditioner to decide the
         % ordering, the rest are to be given to the mech preconditioner
         % for the ordering
         mechBlocks = activeBlockIds(~isFlow);
         flowBlocks = activeBlockIds(isFlow);
         mechOrder = obj.MechPrec.getBlockOrdering(blocks, mechBlocks);
         flowOrder = obj.FlowPrec.getBlockOrdering(blocks, flowBlocks);
         
         % Set the final ordering as mechanics first flow second
         blockOrder = [mechOrder, flowOrder];
      end

      % Set diagonal block scope and propagate to mechanics and flow children
      function setBlockScope(obj, sysInfo, blockIds)
         setBlockScope@preconditioner(obj, sysInfo, blockIds);
         
         % Get all mechanical blocks (also interfaces)
         uBlk = obj.getPhysicsBlocks("displacements");

         % Get all flux blocks (also interfaces)
         flowBlk = obj.getPhysicsBlocks(["pressure", "saturation"]);
         
         % Set the scope
         obj.MechPrec.setBlockScope(sysInfo, obj.precBlocks(uBlk));
         obj.FlowPrec.setBlockScope(sysInfo, obj.precBlocks(flowBlk));
      end

      % Scale the system matrix and RHS to rate form
      function [A, b] = scaleSystem(obj, A, b)
         dt = obj.problemsolver.dt;
         flowRows = obj.getVarBlocks(["pressure", "saturation"]);

         % Scale flow rows of matrix A by 1/dt
         for r = flowRows
            for c = 1:size(A, 2)
               if nnz(A{r, c}) > 0
                  A{r, c} = (1 / dt) * A{r, c};
               end
            end
         end

         % Scale flow rows of RHS b by 1/dt
         for r = flowRows
            b{r} = (1 / dt) * b{r};
         end
      end
   end

   methods (Access = private)
      function isLin = isMechLinear(obj)
         % Set to true, if one domain says it is not then exit
         isLin = true;

         domains = obj.getVarBlocks("displacements");
         interfaces = intersect(obj.getPhysicsBlocks("displacements"),obj.getKindBlocks("interface"));
         
         % if obj.isParentAugmentationNonlinear
         %    return;
         % end

         % Loop over the domains
         for i = interfaces
            % Get the interface solver
            mSolver = obj.problemsolver.interfaces(obj.blockInfo(i).interfId);

            % Ask if it is currently linear
            isLin = mSolver{1}.isLinear();

            % Fast exit if false
            if isLin == false
               return;
            end
         end

         % Loop over the domains
         for d = domains
            % Get the physics solver
            dom = obj.domain(obj.blockInfo(d).domainId);
            mSolver = dom.getPhysicsSolver(obj.blockInfo(d).solverType);
            
            % Ask if it is currently linear
            isLin = mSolver.getMechSolver().isLinear();
            
            % Fast exit if false
            if isLin == false
               return;
            end
         end
      end
   end
end
