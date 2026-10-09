classdef efemPrec < preconditioner

   properties (SetAccess = private,GetAccess=public)

      % Inner preconditioner for mechanics
      innerMech = []

      % A22 block Inverter
      FSAI = []

      % Problemsolver params
      problemsolver

      % Discretizer
      domain
      
      % Number of domains/interfaces
      nDom
      nInt

      % Schur complement matrix on which inner preconditioner was computed on
      S

      % Nonsymmetric matrix tolerance
      nsyTol

      % Regularization factors
      zeroReg = 1e-12
      diagReg = 1e-11
   end

   properties (GetAccess = public,SetAccess = private)
      % Symmetry of the matrix on which the preconditioner has been
      % computed
      PrecSym = true

      % Params structure
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
            error('Not enough arguments for efemPrec apply Left');
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
         
         % Apply inner preconditioner to complete the block upper Gauss Seidel
         x1 = obj.innerMech.ApplyLeft(x11,obj.S);

         % Apply the coupling correction to the state block
         x2 = invC(b2 - B2*x1);
     
         % Compose the solution
         x = [x1;x2];
         
      end

      % Getter for the function handle to apply the right preconditioner
      function x = ApplyRight(obj,b,varargin)
         if nargin < 2
            error('Not enough arguments for efemPrec apply Right');
         end

         x = b;
      end

      % Make A22 slightly more diagonally dominant
      function A22 = regularizeA22(obj,A22)
         diagShift = obj.diagReg * abs(diag(A22));
         diagShift(diagShift == 0) = obj.zeroReg;

         A22 = A22 + diag(diagShift);
      end
      
      function updateStateBlocks(obj, A)
         % Cheap update for reuse

         % Check symmetry
         A22 = A{2,2};
         symm = norm(A22-A22','f')/norm(A22,'f') < obj.nsyTol;

         % Make A22 slightly more diagonally dominant
         A22 = obj.regularizeA22(A22);

         % Compute the factorization of the 22 block
         obj.FSAI.Compute(A22,symm,true);
      
         % Compute reverse Schur complement of A11
         invC = @(x) obj.FSAI.ApplyLeft(x);
         B1 = A{1,2};
         B2 = A{2,1};
      
         % Update function handle while keeping the frozen inner hierarchy S
         obj.Apply_L = @(x) obj.ApplyLeft(x, B1, B2, invC);
      end

      % Constructor Function
      function obj = efemPrec(debugflag,problemsolver,nsyTol,innerPrec)

         % Call the constructor of the abstract class
         obj = obj@preconditioner();
         
         % Set the debugflag
         obj.DEBUGflag = debugflag;

         % Default inner preconditioner to amg if not specified
         if nargin < 4 || isempty(innerPrec)
            innerPrec = 'amg';
         end

         % Get the domains
         obj.problemsolver = problemsolver;
         obj.domain = problemsolver.domains;

         obj.nInt = problemsolver.nInterf;
         obj.nDom = problemsolver.nDom;

         % Create the inner preconditioner for mechanics
         obj.innerMech = createInnerPrec(innerPrec,debugflag,problemsolver,"displacements");

         % Create the inner FSAI solver for fracture jump block
         obj.FSAI = aFSAI(debugflag);

         % Get the parameters
         obj.maxThreads = obj.innerMech.maxThreads;
         if isprop(obj.innerMech, 'params')
            obj.params = obj.innerMech.params;
         end
         obj.nsyTol = nsyTol;
      end
   end
end
