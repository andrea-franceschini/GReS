classdef SAM < preconditioner
%   SAM adaptive Sparse Approximate Mapping.
%
%   This class implements an adaptive Sparse Approximate Mapping (SAM)
%   derived from the abstract 'preconditioner' base class.
%   It is responsible for:
%     - Storing all information necessary for the SAM computation.
%     - Exposing a Compute method to build the SAM for a given
%       couple of matrices and Apply_L / Apply_R handles to apply the SAM + 
%       preconditioner.
%     - If not asked to be used returns the same as if it was not even
%       there.
%     - If obj.N (the SAM matrix) is empty means that is either not used 
%       (useSAM = false) or is the identity matrix (useSAM = true)
%
%   The constructor initializes defaults, enforces thread limits based on
%   system capabilities, and merges user-specified parameters. The class
%   keeps implementation details private and presents a simple interface
%   for preconditioner creation and application.

   properties (GetAccess = public,SetAccess = ?convStrat)
      % Flag to use SAM
      useSAM

      % Preconditioner
      N = []

      % Initial matrix 
      A0
      
      % Symmetry of the matrix on which the preconditioner has been
      % computed
      PrecSym = true

      % Default params for SAM
      nstep = 5
      stp_size = 1
      epss = 1e-4
      
      % Max Threads
      maxThreads

      % Recomputing
      alpha = 0.5
      requestComp = false % starts at false and modified inside convStrat
      precJustComputed = true
      firstCompPercDegrad = 0.01
      nSolveSinceLastComp
      firstSolveTAfterComp = 0

      % Statistics
      aTimeComp = 0
      nComp = 0
      Delta_T = zeros(2,1)
      cumTSolveAfterComp = 0
      tSetup

   end
   properties (Access = public)
      CompLin = []
   end

   methods (Access = public)

      % Function to compute the preconditioner
      Compute(obj,A,sym,varargin)

      % Getter for the function handle to apply the left preconditioner
      function x = ApplyLeft(obj,b,varargin)
         if nargin < 3
            error('Not enough arguments for SAM apply Left');
         end
         Prec = varargin{1};

         % Sanity check
         if isempty(obj.N)
            error('SAM needs to be non empty');
         end
         
         x = sam_apply_left(obj.N, Prec.Apply_L, b);
      end

      % Getter for the function handle to apply the right preconditioner
      function x = ApplyRight(obj,b,varargin)
         if nargin < 3
            error('Not enough arguments for SAM apply Right');
         end
         Prec = varargin{1};

         x = Prec.Apply_R(b);
      end

      % Constructor Function
      function obj = SAM(useSAM)

         % Call the constructor of the abstract class
         obj = obj@preconditioner();

         % Set if to use the SAM
         obj.useSAM = useSAM;

         % Set maximum number of threads to use if the system provides less
         obj.maxThreads = maxNumCompThreads;

      end

      % Reset after preconditioner computation
      function reset(obj)
         % Reset only if needed
         if obj.useSAM
            obj.N = [];
            obj.requestComp = false;
            obj.nSolveSinceLastComp = 0;
            obj.precJustComputed = true;
         end
      end

      % Get old matrix
      function getOldMat(obj,A0,symm)
         % Save only if needed
         if obj.useSAM
            obj.A0 = A0;
            obj.PrecSym = symm;
         end
      end
   end
end


