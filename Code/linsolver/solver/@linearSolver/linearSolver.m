classdef linearSolver < handle
% linearSolver - Handle class that manages linear solver selection, preconditioner
%                creation (Chronos_Lab), solving strategy and statistics.
%
% Usage:
%   obj = linearSolver(generalsolver, physname)
%
% Inputs:
%   generalsolver - struct or object containing simulation parameters and
%                   interfaces required to configure solver and preconditioner.
%   physname      - string with the physics name used to select appropriate
%                   preconditioner setup.
%
% Description:
%   This class detects presence of the Chronos_Lab third-party library and,
%   when available and compiled, configures a Chronos preconditioner and a
%   solver strategy (GMRES, etc.). It also supports a MATLAB fallback if
%   Chronos is missing or not compiled. The object gathers statistics about
%   preconditioner computation, solve times and iteration counts, and can
%   print a summary via printStats.
%
% Main methods:
%   linearSolver(generalsolver, physname) - constructor: checks Chronos,
%       sets up preconditioner and solver defaults.
%   printStats(varargin) - print collected statistics; 'short' option prints
%       only summary.
%
% Notes:
%   - If generalsolver.simparams.linSolverParams.useMatlab is set, constructor
%     will exit early to allow MATLAB solvers only.
% Brief property descriptions:
% DEBUGflag       - boolean: enable debug logging and extra checks.
% matlabMaxSize   - numeric: maximum matrix size for which MATLAB fallback is enforced.
% nsyTol          - numeric: tolerance used to detect near-symmetric matrices.
% convStrat       - object: convergence strategy handler (convStrat class instance).
% ChronosFlag     - boolean: true if Chronos_Lab preconditioner support is available.
% requestPrecComp - boolean: flag requesting that the preconditioner be (re)computed.
% alpha           - numeric: scaling parameter used during recomputation of the preconditioner.
% x0              - vector: starting vector for iterative solvers.
% SolverType      - string: name of the iterative solver in use (e.g. 'gmres').
% Prec            - object: preconditioner instance (Chronos preconditioner wrapper).
% iterConfigOld   - numeric: stored iteration configuration version to detect changes.
% whenComputed    - vector: physical times at which the preconditioner was computed.
% aTimeComp       - numeric: accumulated time spent computing preconditioners.
% aTimeSolve      - numeric: accumulated time spent in solver calls.
% nSolve          - integer: number of times the linear solver was invoked.
% nComp           - integer: number of preconditioner computations performed.
% maxIter         - integer: maximum iterations observed across solves.
% aIter           - numeric: accumulated iterations (for averaging/stats).
% cumTSolveAfterPrec - numeric: cumulative solve time measured after preconditioner build.
% iterLin         - vector: recorded iteration counts per linear solve.
% solveTLin       - vector: recorded solve times per linear solve.
% symFlagLin      - vector: flags indicating symmetry detection per solve.
% precCompLin     - vector: flags indicating whether prec was computed for each solve.
% newtonLin       - vector: Newton iteration indices associated with each linear solve.
% timeLin         - vector: physical times corresponding to each linear solve.
% Delta_T         - vector: time step sizes used when recording solves.
% useSAM          - boolean: whether SAM monitoring/analysis is enabled.
% SAM             - object: SAM instance.
% params          - struct: runtime parameters for the solver (tolerances, maxit, restart, etc.)

   properties (SetAccess = ?convStrat, GetAccess = public)

      % Flag for debug
      DEBUGflag = false
      matlabMaxSize = 2e4

      % Utils flags
      nsyTol = 100*eps

      % Convergence strategy handler
      convStrat
      
      % Flag for Chronos existance
      ChronosFlag = false

      % Flag to request Preconditioner computation
      requestPrecComp = true
      alpha = 1

      % starting vector
      x0 = []

      % Solver Type
      SolverType

      % Preconditioner object
      Prec

      % General solver
      generalsolver
      iterConfigOld = 1

      % Statistics
      whenComputed = []
      aTimeComp = 0
      aTimeSolve = 0
      nSolve = 0
      nComp = 0
      maxIter = -1
      aIter = 0
      cumTSolveAfterPrec = 0
      iterLin = []
      solveTLin = []
      symFlagLin = []
      precCompLin = []
      newtonLin = []
      timeLin = []
      Delta_T = []

      % SAM object
      useSAM = false
      SAM = []

      % Params struct
      params
   end

   methods (Access = public)

      % Constructor Function
      function obj = linearSolver(generalsolver,physname)

         % Check if chronos is available
         ChronosDir = fullfile(gres_root,'ThirdPartyLibs','Chronos_Lab','sources');

         % Possible the user wants to use matlab even if the size is sufficient
         if isfield(generalsolver.simparams.linSolverParams, 'useMatlab')
           if generalsolver.simparams.linSolverParams.useMatlab == 1
             gresLog().log(3,'The user requested the use of matlab\n');
             return;
           end
         end

         if isfolder(ChronosDir)
           if strcmp(computer('arch'),'maca64')
             fileMex = fullfile(ChronosDir,'Preconditioner','AMG','filter','MEX_Prol_Filter','FilterProl_wrap.mexmaca64');
           elseif strcmp(computer('arch'),'win64')
             fileMex = fullfile(ChronosDir,'Preconditioner','AMG','filter','MEX_Prol_Filter','FilterProl_wrap.mexw64');
           else
             fileMex = fullfile(ChronosDir,'Preconditioner','AMG','filter','MEX_Prol_Filter','FilterProl_wrap.mexa64');
           end

            if ~isfile(fileMex)
               gresLog().warning(1,'Chronos_Lab submodule is present, but not compiled. Using matlab fallback');
               return;
            end

            obj.generalsolver = generalsolver;

            % Create the preconditioner object, check if the physics is supported
            [obj.Prec,obj.ChronosFlag] = obj.choosePrec(obj.DEBUGflag,generalsolver,physname);

            % Create SAM object, deals with everything inside if not asked
            % to be used
            obj.SAM = SAM(obj.useSAM);

            % Non supported physics for the preconditioner
            if ~obj.ChronosFlag
               return;
            end

            % Chronos exists
            addpath(genpath(ChronosDir));

            % First time solving request preconditioner computation
            obj.params.iter = -1;
            obj.params.lastRelres = 1e10;

            % Choose the relative tolerance strategy
            obj.convStrat = convStrat(generalsolver);

            % Get default values
            chronos_xml_default = fullfile(gres_root,'Code','linsolver','XML_setup','chronos_xml_setup.xml');

            % Read Defaults
            data = readInput(chronos_xml_default);

            % Get the solver type
            obj.SolverType = lower(data.solver);
            obj.params.maxit = 500;

            % if GMRES get restart value
            if strcmp(obj.SolverType,'gmres')
               obj.params.restart = data.general.restart;
            else
               obj.params.restart = 100;
            end
         end
      end

      function printStats(obj,varargin)
         % Allocate for stats if sam not used
         if ~obj.useSAM
            obj.SAM.CompLin = zeros(size(obj.precCompLin));
         end
         
         % Check if the user wants only the summary review
         if ~isempty(varargin)
            string = varargin{1};
         else
            string = [];
         end

         % The user did not ask for the summary only, show full info
         if ~strcmpi(string,'short')
            % Print the time steps at which the preconditioner was computed
            fprintf('The preconditioner was computed at time(s):\n');
            for i = 1:length(obj.whenComputed)
               fprintf('             %d\t%e\n',i,obj.whenComputed(i));
            end

            fprintf('\n----------------------------------------------------------------------\n')
            fprintf('| %8s | %6s | %4s | %8s | %7s | %8s | %8s |\n','PhysTime','Sol N.','Iter','SolTime','Symm','PrecTime','DeltaT');
            fprintf('----------------------------------------------------------------------\n')
            
            for i = 1:size(obj.solveTLin,2)
               fprintf('| %.2e | %6d | %4d | %.2e | %.1e | %.2e | %.2e |\n',obj.timeLin(i),i,obj.iterLin(i),obj.solveTLin(i),obj.symFlagLin(i),obj.precCompLin(i)+obj.SAM.CompLin(i),obj.Delta_T(i));
            end
         end

         fprintf('Average Preconditioner computation time = %e\n',(obj.aTimeComp/obj.nComp));
         fprintf('Preconditioner computed %d times\n',length(obj.whenComputed));
         if obj.useSAM
            nnzSAMcomp = size(nonzeros(obj.SAM.CompLin),1);
            fprintf('Average SAM computation time = %e\n',sum(obj.SAM.CompLin)/nnzSAMcomp);
            fprintf('SAM computed %d times\n',nnzSAMcomp);
         end
         fprintf('Average Solve time = %e\n',(obj.aTimeSolve/obj.nSolve));
         fprintf('Average number of iterations = %f\n',(obj.aIter/obj.nSolve));
         fprintf('Max number of iterations = %d\n',obj.maxIter);
         fprintf('Total time for computation of the linear systems = %e\n',obj.aTimeComp+obj.aTimeSolve+sum(obj.SAM.CompLin));
         % fprintf('Used %d threads during mex\n',obj.Prec.maxThreads);

      end

      % Function to solve the system
      [x,flag] = SolveLin(obj,A,b,time,nonlinIter,isLinear);

      % Function to get if chronos is to be used and if so instanciate the
      % preconditioner
      [Prec,ChronosFlag] = choosePrec(obj,debugflag,problemsolver,physname);

      % Function to get the total time taken by the linear solver for
      % preconditioner computation and solve step
      function tot = getTotalTime(obj)
         if ~obj.useSAM || ~obj.ChronosFlag
            obj.SAM.CompLin = zeros(size(obj.precCompLin));
         end
         tot = obj.aTimeComp+obj.aTimeSolve+sum(obj.SAM.CompLin);
      end

   end
end
