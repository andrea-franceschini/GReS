classdef aFSAI < preconditioner
%   aFSAI adaptive Factorized Sparse Approximate Inverse.
%
%   This class implements an adaptive Factorized Sparse Approximate Inverse (aFSAI)
%   preconditioner derived from the abstract 'preconditioner' base class.
%   It is responsible for:
%     - Constructing and storing the internal preconditioner operator(s),
%       along with any symmetry information and threading limits.
%     - Exposing a Compute method to build the preconditioner for a given
%       matrix and Apply_L / Apply_R handles to apply the preconditioner.
%
%   The constructor initializes defaults, enforces thread limits based on
%   system capabilities, and merges user-specified parameters. The class
%   keeps implementation details private and presents a simple interface
%   for preconditioner creation and application.

   properties (GetAccess = public,SetAccess = private)
       % Preconditioner
      Prec = []
      
      % Symmetry of the matrix on which the preconditioner has been
      % computed
      PrecSym = true
      
      % Params
      param
      maxThreads
      nstep = 10
      step_size = 1
      epsilon = 1e-5

   end

   methods (Access = public)

      % Function to compute the preconditioner
      Compute(obj,A,sym,varargin)

      % Getter for the function handle to apply the left preconditioner
      function x = ApplyLeft(obj,b,varargin)
         if nargin < 2
            error('Not enough arguments for FSAI apply Left');
         end

         % FSAI application
         x = obj.Prec.omega*obj.Prec.right*(obj.Prec.left*b);  
      end

      % Getter for the function handle to apply the right preconditioner
      function x = ApplyRight(obj,b,varargin)
         if nargin < 2
            error('Not enough arguments for FSAI apply Right');
         end

         x = b;
      end
         

      % Constructor Function
      function obj = aFSAI(debugflag)

         % Call the constructor of the abstract class
         obj = obj@preconditioner();

         % Set the debugflag
         obj.DEBUGflag = debugflag;

         % Set maximum number of threads to use if the system provides less
         obj.maxThreads = maxNumCompThreads;

         % Get the different parameters 
         obj.param.nthread      = obj.maxThreads;
         obj.param.nstep        = obj.nstep;
         obj.param.step_size    = obj.step_size;
         obj.param.epsilon      = obj.epsilon;
         obj.param.method       = 'afsai_sym';
      end
   end
end


