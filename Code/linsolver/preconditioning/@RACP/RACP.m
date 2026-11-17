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
%     - Provides Compute method for building the operator
%       
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

   properties (SetAccess = private, GetAccess = public)
      % Inner preconditioner for augmented A11 block (e.g. aAMG, fixedStress)
      innerPrec
   end

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
      A = Compute(obj,A,symMat,varargin)

      % Set diagonal block scope and propagate variable blocks to inner preconditioner
      function setBlockScope(obj, sysInfo, blockIds)
         setBlockScope@preconditioner(obj, sysInfo, blockIds);
         obj.innerPrec.setBlockScope(sysInfo, obj.precBlocks(obj.getKindBlocks("variable")));
      end

      % Set inner sub-preconditioner
      function setSubPrecs(obj, varargin)
         if isa(varargin{1}, 'preconditioner')
            obj.innerPrec = varargin{1};
         else
            obj.innerPrec = varargin{2};
         end
      end

      % Determine adaptive block ordering for activeBlockIds
      function blockOrder = getBlockOrdering(obj, blocks, activeBlockIds)
         activeMeta = blocks(activeBlockIds);
         isSaddle = [activeMeta.kind] == "interface";
         mainBlocks = activeBlockIds(~isSaddle);
         saddleBlocks = activeBlockIds(isSaddle);
         mainOrder = obj.innerPrec.getBlockOrdering(blocks, mainBlocks);
         blockOrder = [mainOrder, saddleBlocks];
      end

      % Getter for the function handle to apply the left preconditioner
      function x = ApplyLeft(obj,b,varargin)
         if numel(varargin) >= 3
            A12     = varargin{1};
            A21     = varargin{2};
            inv_D22 = varargin{3};
         else
            error('Not enough arguments for RACP apply Left');
         end
         
         n1 = size(A12, 1);
         b1 = b(1:n1);
         b2 = b(n1+1:end);

         % Augment RHS for A11 block
         rhs1 = b1 + A12 * (inv_D22 * b2);

         % Solve augmented block polymorphically via innerPrec.Apply_L
         x1 = obj.innerPrec.Apply_L(rhs1);

         % Backsubstitute for interface multipliers
         x2 = inv_D22 * (b2 - A21 * x1);

         x = [x1; x2];
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
         obj.innerPrec = aAMG(debugflag,problemsolver,physname);

         obj.maxThreads = obj.innerPrec.maxThreads;
         obj.params = obj.innerPrec.params;
         
      end

      % Compute local Augmented matrix
      [A11_aug,inv_D22,A12,A21] = cpt_localAug(obj,A,symm,varargin)
   
      % Condenses the domains and interfaces in a 2x2 matrix
      A = condenseDomains(obj,A)
   end
end
