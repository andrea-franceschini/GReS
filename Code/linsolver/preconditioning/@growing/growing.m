classdef growing < preconditioner
%   Growing Preconditioner
%
%   This class implements the growing as a subclass of the abstract preconditioner
%   class. This preconditioner leverages on an internal AMG preconditioner
%   which is progressively grown using a jacobi diagonal preconditioner.
%   Once needed the full AMG is then recomputed on the whole new domain
%
%   Key features:
%     - Holds an inner AMG preconditioner (AMG) 
%
%   Usage:
%     obj = growing(debugflag,problemsolver,physname)
%       Constructs the growing preconditioner with debugging control,
%       a handle to the problemsolver (which must expose domain info),
%       and a string describing the physics (e.g., 'pressure',
%       'displacements', 'displacements_contact').
%
%   Notes:
%     - The class expects the problemsolver.domains to be available and
%       that the provided physname matches supported physics types.
%     - Inner AMG behavior and parameters are encapsulated in the aAMG
%       instance created during construction.

   properties (GetAccess = public,SetAccess = private)

      % Inner preconditioners
      AMG

      % Problemsolver params
      problemsolver

      % Discretizer
      domain

      % Symmetry of the matrix on which the preconditioner has been
      % computed
      PrecSym = true

      % AMG params structure
      params

      % Physics
      phys

      % Max Threads
      maxThreads

      % Size difference from when the preconditioner was originally
      % computed
      sizeDiff = 0
      
   end

   methods
      % Function to check if the system has grown
      x0 = checkGrowth(obj, linsolver, b);

      % Function to compute the growing preconditioner
      function Compute(obj,A,symMat,varargin)
      
         % Set the symmetry of the full preconditioner
         obj.PrecSym = symMat;

         % Reset size diff to 0 after full recomputation
         obj.sizeDiff = 0;

         % Check that the matrix is in sparse format and not in cell format
         if iscell(A)
            A = A{1,1};
         end

         % Compute the test space
         TV0 = ones(size(A,1),1);
      
         % Compute the amg for block 11
         obj.AMG.Compute(A,obj.PrecSym,TV0,false);
      
         obj.Apply_L = @(x) obj.ApplyLeft(x,A);
         obj.Apply_R = @(x) obj.ApplyRight(x);
      end

      % Update the application of the preconditioner if the size has
      % changed
      function updateGrowingPrec(obj,Amat)
      
         % Check if the size is different from when the preconditioner was
         % computed
         if obj.sizeDiff > 0 
            invD = 1./diag(Amat(end-obj.sizeDiff+1:end,end-obj.sizeDiff+1:end));
            obj.Apply_L = @(x) [obj.AMG.Apply_L(x(1:end-obj.sizeDiff)); 
                                invD.*x(end-obj.sizeDiff+1:end)];
         end
      end

      % Getter for the function handle to apply the left preconditioner
      function x = ApplyLeft(obj,b,varargin)
         if nargin < 3
            error('Not enough arguments for growing apply Left');
         end

         A = varargin{1};

         x = obj.AMG.ApplyLeft(b,A);
      end

      % Getter for the function handle to apply the right preconditioner
      function x = ApplyRight(obj,b,varargin)
         if nargin < 2
            error('Not enough arguments for growing apply Right');
         end

         x = b;
      end

      % Constructor Function
      function obj = growing(debugflag,problemsolver,physname)

         % Call the constructor of the abstract class
         obj = obj@preconditioner();
         
         % Set the debugflag
         obj.DEBUGflag = debugflag;

         % Get the domains
         obj.problemsolver = problemsolver;
         obj.domain = problemsolver.domains;

         if problemsolver.nInterf ~= 0
            error("Growing preconditioner does not support interfaces")
         end
         
         if problemsolver.nDom ~= 1
            error("Growing preconditioner does not support multiple domains")
         end

         % Supported Single Physics
         if(contains(physname,'pressure') || contains(physname,'u'))
            obj.phys = 0;
         else
            disp(physname);
            error('Non supported Physics for growing preconditioner');
         end

         % Create the inner AMG preconditioner
         obj.AMG = aAMG(debugflag,problemsolver,physname);

         obj.maxThreads = obj.AMG.maxThreads;
         obj.params = obj.AMG.params;
      end
   end
end