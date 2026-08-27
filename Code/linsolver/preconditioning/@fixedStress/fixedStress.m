classdef fixedStress < preconditioner

   properties (Access = private)

      % Inner preconditioners
      AMGMech = []
      AMGFlux = []

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
      Compute(obj,A,sym,varargin)

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

         % Apply AMG for fluid part
         x2 = obj.AMGFlux.ApplyLeft(b(n1+1:end),S);

         % Apply the top part of the block preconditioner
         x11 = b(1:n1) - B1*x2;
         x1 = obj.AMGMech.ApplyLeft(x11,A11);

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
      function obj = fixedStress(debugflag,problemsolver)

         % Call the constructor of the abstract class
         obj = obj@preconditioner();
         
         % Set the debugflag
         obj.DEBUGflag = debugflag;

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

         % Create the inner AMG preconditioners
         obj.AMGFlux = aAMG(debugflag,problemsolver,"pressure");
         obj.AMGMech = aAMG(debugflag,problemsolver,"displacements");

         obj.maxThreads = obj.AMGFlux.maxThreads;
         obj.params = obj.AMGFlux.params;
      end

      % Function for treating the dirichlet boundary conditions
      A = treatDirBC(obj,A,sym)

   end
end