classdef RACP < preconditioner
%   RACP Reverse Augmented Constraint Preconditioner
%
%   This class implements the RACP (Reverse Augmented Constraint
%   Preconditioner) as a subclass of the abstract preconditioner
%   class. RACP provides a block preconditioning strategy for coupled
%    problems arising from finite element discretizations. 
%    It builds a reverse-augmented constraint framework while delegating 
%   inner solves to an algebraic multigrid (AMG) preconditioner.
%
%   Key features:
%     - Holds an inner AMG preconditioner (AMG) for approximate inversion
%       of block diagonal or Schur-complement approximations.
%     - Stores references to the problem solver and domain discretization,
%       enabling assembly and treatment of domain/interface quantities.
%     - Supports different physics modes (pressure, displacement, contact)
%       and adapts preconditioning strategy accordingly.
%     - Provides Compute and treatDirBC methods for building the operator
%       and handling Dirichlet boundary conditions.
%
%   Usage:
%     obj = RACP(debugflag,problemsolver,physname)
%       Constructs the RACP preconditioner with debugging control,
%       a handle to the problemsolver (which must expose domain info),
%       and a string describing the physics (e.g., 'pressure',
%       'displacements', 'displacements_contact').
%
%   Notes:
%     - The class expects the problemsolver.domains to be available and
%       that the provided physname matches supported physics types.
%     - Inner AMG behavior and parameters are encapsulated in the aAMG
%       instance created during construction.

   properties (Access = private)

      % Inner preconditioners
      AMG

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
      Compute(obj,A,symMat,varargin)

      % % Function to apply the multidomain RACP procedure
      % y = applyMultiRACP(obj,A11_aug,B1,B2,inv_D22,x);

      % Getter for the function handle to apply the left preconditioner
      function x = ApplyLeft(obj,b,varargin)
         if nargin < 6
            error('Not enough arguments for RACP apply Left');
         end
         A11_aug = varargin{1};
         A12     = varargin{2};
         A21     = varargin{3};
         inv_D22 = varargin{4};
         
         x = apply_RevAug(obj.AMG.Prec,A11_aug,A12,A21,inv_D22,b);
      end

      % Getter for the function handle to apply the right preconditioner
      function x = ApplyRight(obj,b,varargin)
         if nargin < 2
            error('Not enough arguments for RACP apply Right');
         end

         x = b;
      end

      % Constructor Function
      function obj = RACP(debugflag,problemsolver,physname)

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

         % Create the inner AMG preconditioner
         obj.AMG = aAMG(debugflag,problemsolver,physname);

         obj.maxThreads = obj.AMG.maxThreads;
         obj.params = obj.AMG.params;
         
      end

      % Function for treating the dirichlet boundary conditions
      A = treatDirBC(obj,A,sym)

      % Compute local Augmented matrix for multidom == false
      [A11_aug,inv_D22] = cpt_localAug(obj,A11,A12,A21,A22,symm)

   end
end