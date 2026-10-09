classdef fixedStress < preconditioner

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
      % Inner preconditioners
      innerMech = []
      innerFlux = []

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
            error('Not enough arguments for fixedStress apply Left');
         end

         S = varargin{1};
         A11 = varargin{2};
         B1 = varargin{3};

         % Get mechanics size
         n1 = size(A11,1);

         % Apply inner preconditioner for fluid part
         x2 = obj.innerFlux.ApplyLeft(b(n1+1:end),S);

         % Apply the top part of the block preconditioner
         x11 = b(1:n1) - B1*x2;
         x1 = obj.innerMech.ApplyLeft(x11,A11);

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

      % Constructor Function
      function obj = fixedStress(debugflag,problemsolver,innerPrec,innerPrecFlux)

         % Call the constructor of the abstract class
         obj = obj@preconditioner();
         
         % Set the debugflag
         obj.DEBUGflag = debugflag;

         % Default inner preconditioners to amg if not specified
         if nargin < 3 || isempty(innerPrec)
            innerPrec = 'amg';
         end
         if nargin < 4 || isempty(innerPrecFlux)
            innerPrecFlux = innerPrec;
         end

         % Get the domains
         obj.problemsolver = problemsolver;
         obj.domain = problemsolver.domains;

         obj.nInt = problemsolver.nInterf;
         obj.nDom = problemsolver.nDom;

         % Check the number of interfaces and domains
         if obj.nInt ~= 0
            interfacein = problemsolver.interfaces;
         else
            interfacein = {};
         end

         % Create the inner preconditioners
         obj.innerFlux = createInnerPrec(innerPrecFlux,debugflag,problemsolver,"pressure");
         obj.innerMech = createInnerPrec(innerPrec,debugflag,problemsolver,"displacements");

         obj.maxThreads = obj.innerFlux.maxThreads;
         if isprop(obj.innerFlux, 'params')
            obj.params = obj.innerFlux.params;
         end
      end
   end
end
