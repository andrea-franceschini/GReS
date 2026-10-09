classdef RACP < preconditioner
%   RACP Reverse Augmented Constraint Preconditioner
%
%   This class implements the RACP (Reverse Augmented Constraint
%   Preconditioner) as a subclass of the abstract preconditioner
%   class. RACP provides a block preconditioning strategy for coupled
%   problems arising from finite element discretizations. 
%   It builds a reverse-augmented constraint framework while delegating 
%   inner solves to an inner preconditioner (e.g. AMG or FSAI).
%
%   Key features:
%     - Holds an inner preconditioner (innerPrec) for approximate inversion
%       of block diagonal or Schur-complement approximations.
%     - Stores references to the problem solver and domain discretization,
%       enabling assembly and treatment of domain/interface quantities.
%     - Supports different physics modes (pressure, displacement, contact)
%       and adapts preconditioning strategy accordingly.
%     - Provides Compute method for building the operator
%
%   Usage:
%     obj = RACP(debugflag,problemsolver,physname,innerPrec)
%       Constructs the RACP preconditioner with debugging control,
%       a handle to the problemsolver (which must expose domain info),
%       a string describing the physics (e.g., 'pressure',
%       'displacements', 'displacements_contact'), and an optional
%       inner preconditioner type ('amg' or 'fsai').

   properties (Access = private)

      % Problemsolver params
      problemsolver

      % Discretizer
      domain
      
      % Number of domains/interfaces
      nDom
      nInt

      % RACP Gamma
      gamma = 1.0

   end

   properties (GetAccess = public,SetAccess = private)
      % Inner preconditioner
      innerPrec

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
      A = Compute(obj,A,symMat,varargin)

      % Getter for the function handle to apply the left preconditioner
      function x = ApplyLeft(obj,b,varargin)
         if nargin < 6
            error('Not enough arguments for RACP apply Left');
         end
         A11_aug = varargin{1};
         A12     = varargin{2};
         A21     = varargin{3};
         inv_D22 = varargin{4};
         
         % Get domain augmented block size
         n1 = size(A11_aug,1);

         % Partition the right-hand side vector
         x1 = b(1:n1,:);
         x2 = b(n1+1:end,:);

         % Compute augmented rhs for domain block
         b1 = x1 + A12*(inv_D22*x2);

         % Apply inner preconditioner to domain augmented block
         y1 = obj.innerPrec.ApplyLeft(b1,A11_aug);

         % Compute interface block solution
         y2 = inv_D22*(A21*y1 - x2);

         % Compose the solution
         x = [y1; y2];
      end

      % Getter for the function handle to apply the right preconditioner
      function x = ApplyRight(obj,b,varargin)
         if nargin < 2
            error('Not enough arguments for RACP apply Right');
         end

         x = b;
      end

      % Constructor Function
      function obj = RACP(debugflag,problemsolver,physname,innerPrec)

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

         % Check the number of interfaces and domains
         if obj.nInt ~= 0
            interfacein = problemsolver.interfaces;
         else
            interfacein = {};
         end

         % Supported Single Physics
         if(contains(physname,'pressure') || contains(physname,'u'))
            obj.phys = 0;
         elseif(contains(physname,'displacements')) 
            obj.phys = 1;
            % Check if there is contact 
            if any(cellfun(@(o) isa(o,'SolidMechanicsContact'),interfacein))
               obj.phys = 1.1;
            end
         else
            disp(physname);
            error('Non supported Physics for preconditioner');
         end

         % Create the inner preconditioner
         obj.innerPrec = createInnerPrec(innerPrec,debugflag,problemsolver,physname);

         obj.maxThreads = obj.innerPrec.maxThreads;
         if isprop(obj.innerPrec, 'params')
            obj.params = obj.innerPrec.params;
         end
         
      end

      % Compute local Augmented matrix
      [A11_aug,inv_D22] = cpt_localAug(obj,A11,A12,A21,A22,symm)
   
      % Condenses the domains and interfaces in a 2x2 matrix
      A = condenseDomains(obj,A)
   end
end
