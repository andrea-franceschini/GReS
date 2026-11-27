classdef efemPrec < preconditioner

   properties (SetAccess = private,GetAccess=public)

      % Inner preconditioner
      innerPrec = []

      % A22 block Inverter
      FSAI = []

      % Problemsolver params
      problemsolver

      % Discretizer
      domain
      
      % Number of domains/interfaces
      nDom

      % Schur complement matrix on which AMG was computed on
      S

      % Nonsymmetric matrix tolerance
      nsyTol
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
         x1 = obj.innerPrec.ApplyLeft(x11,obj.S);

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

      function updateStateBlocks(obj, A)
         % Cheap update for reuse

         uBlk = obj.getVarBlocks("displacements");
         jBlk = obj.getVarBlocks("fractureJump");

         % Check symmetry
         A22 = A{jBlk,jBlk};
         symm = norm(A22-A22','f')/norm(A22,'f') < obj.nsyTol;

         % Compute the factorization of the 22 block
         obj.FSAI.Compute(A22,symm,true);
         
         % Compute reverse Schur complement of A11
         invC = @(x) obj.FSAI.ApplyLeft(x);
         B1 = A{uBlk,jBlk};
         B2 = A{jBlk,uBlk};
      
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
         obj.nDom = problemsolver.nDom;

         % Create the inner AMG preconditioner
         obj.innerPrec = aAMG(debugflag,problemsolver,"displacements");

         % Create the inner FSAI solver
         obj.FSAI = aFSAI(debugflag);

         % Get the parameters
         obj.maxThreads = obj.innerPrec.maxThreads;
         obj.params = obj.innerPrec.params;
         obj.nsyTol = nsyTol;
         obj.linsolverTol = problemsolver.simparams.relTol;
      end


      % Set sub-preconditioners for EFEM (mechanics)
      function setSubPrecs(obj, varargin)
         for k = 1:2:numel(varargin)
            if strcmpi(string(varargin{k}), "displacements")
               obj.innerPrec = varargin{k+1};
            end
         end
      end

      % Determine adaptive block ordering for activeBlockIds
      function blockOrder = getBlockOrdering(obj, blocks, activeBlockIds)
         activeMeta = blocks(activeBlockIds);
         uBlocks = activeBlockIds([activeMeta.varName] == "displacements");
         jBlocks = activeBlockIds([activeMeta.varName] == "fractureJump");
         uOrder = obj.innerPrec.getBlockOrdering(blocks, uBlocks);
         blockOrder = [uOrder, jBlocks];
      end

      % Set diagonal block scope and propagate to AMG and FSAI children
      function setBlockScope(obj, sysInfo, blockIds)
         setBlockScope@preconditioner(obj, sysInfo, blockIds);
         uBlk = obj.getVarBlocks("displacements");
         jBlk = obj.getVarBlocks("fractureJump");
         obj.innerPrec.setBlockScope(sysInfo, obj.precBlocks(uBlk));
         obj.FSAI.setBlockScope(sysInfo, obj.precBlocks(jBlk));
      end

   end
end
