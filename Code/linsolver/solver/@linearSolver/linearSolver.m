classdef linearSolver < handle
% LINEARSOLVER  Wrapper around external (Chronos) and MATLAB linear solvers.
%
%   This class manages optional preconditioning, solver configuration,
%   timing and statistics for preconditioner construction and linear solves.
%
%   PROPERTIES (public get, private set):
%     DEBUGflag        - Flag for debug output (logical)
%     matlabMaxSize    - Threshold size to force MATLAB solver (numeric)
%     nsyTol           - Numerical symmetry tolerance (numeric)
%     ChronosFlag      - True if Chronos preconditioner is available (logical)
%     requestPrecComp  - Request preconditioner computation (logical)
%     x0               - Starting vector for iterative solvers (numeric vector)
%     SolverType       - Solver type string (e.g., 'gmres')
%     Prec             - Preconditioner object (Chronos wrapper)
%     generalsolver    - Reference to the nonlinear solver object
%     iterConfigOld    - Configuration flag for iterative solver reuse
%     whenComputed     - Times at which preconditioner was computed
%     aTimeComp        - Accumulated time spent computing preconditioner
%     aTimeSolve       - Accumulated time spent solving linear systems
%     nSolve           - Number of solves performed
%     nComp            - Number of preconditioner computations performed
%     maxIter          - Maximum iterations observed across solves
%     aIter            - Accumulated iteration counts
%     iterLin          - Per-solve iteration counts
%     solveTLin        - Per-solve solve times
%     symFlagLin       - Per-solve symmetry flags
%     precCompLin      - Per-solve preconditioner computation flags
%     newtonLin        - Per-solve Newton step indices
%     timeLin          - Per-solve timestamps
%     params           - Struct of solver parameters (tol, maxit, restart, ...)
%
%   METHODS:
%     linearSolver(generalsolver, physname)
%       Constructor. Checks for Chronos library and compiled mex, creates
%       the preconditioner object when supported, reads default Chronos XML
%       settings, and initializes solver parameters. Sets ChronosFlag=false
%       if Chronos is unavailable or useMatlab is forced via
%       generalsolver.simparams.linSolverParams.useMatlab.
%
%     printStats()
%       Prints accumulated statistics for preconditioner construction and
%       linear solves. Also prints a per-solve table of timing, iterations, 
%       symmetry flag and preconditioner computation time.
%
%     Solve(obj, A, b, time)
%       Solves the linear system A*x = b using the configured solver and
%       preconditioner. Returns solution x and a convergence flag.
%

   properties (SetAccess = ?convStrat, GetAccess = public)

      % Flag for debug
      DEBUGflag = false
      matlabMaxSize = 2e5

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

      % Ruiz params
      nIterRuiz = 0
      tolRuiz = 1e-3

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
            obj.params.maxit = data.general.maxit;

            % if GMRES get restart value
            if strcmp(obj.SolverType,'gmres')
               obj.params.restart = data.general.restart;
            else
               obj.params.restart = 100;
            end

            % Get the values if ruiz symmetric scaling is asked for by the user
            if isfield(generalsolver.simparams.linSolverParams, 'ruizIter')
               obj.nIterRuiz = generalsolver.simparams.linSolverParams.ruizIter;
            end
            if isfield(generalsolver.simparams.linSolverParams, 'ruizTol')
               obj.tolRuiz = generalsolver.simparams.linSolverParams.ruizTol;
            end
         end
      end

      function printStats(obj,varargin)
         fprintf('\n\n\nAverage Preconditioner computation time = %e\n',(obj.aTimeComp/obj.nComp));
         fprintf('The preconditioner was computed at time(s):\n');
         for i = 1:length(obj.whenComputed)
            fprintf('             %d\t%e\n',i,obj.whenComputed(i));
         end
         fprintf('\nAverage Solve time = %e\n',(obj.aTimeSolve/obj.nSolve));
         fprintf('\nAverage number of iterations = %f\n',(obj.aIter/obj.nSolve));
         fprintf('Max number of iterations = %d\n',obj.maxIter);
         fprintf('\nTotal time for computation of the linear systems = %e\n',obj.aTimeComp+obj.aTimeSolve);

         if ~isempty(varargin)
            string = varargin{1};
         else
            string = [];
         end

         if ~strcmpi(string,'short')
            fprintf('\nUsed %d threads during mex\n',obj.Prec.maxThreads);
            fprintf('\n----------------------------------------------------------------------\n')
            fprintf('| %8s | %6s | %4s | %8s | %7s | %8s | %8s |\n','PhysTime','Sol N.','Iter','SolTime','Symm','PrecTime','DeltaT');
            fprintf('----------------------------------------------------------------------\n')
            
            for i = 1:size(obj.solveTLin,2)
               fprintf('| %.2e | %6d | %4d | %.2e | %.1e | %.2e | %.2e |\n',obj.timeLin(i),i,obj.iterLin(i),obj.solveTLin(i),obj.symFlagLin(i),obj.precCompLin(i),obj.Delta_T(i));
            end
         end
      end

      % Function to solve the system
      [x,flag] = SolveLin(obj,A,b,time,nonlinIter,isLinear);

      % Function to get if chronos is to be used and if so instanciate the
      % preconditioner
      [Prec,ChronosFlag] = choosePrec(obj,debugflag,problemsolver,physname);

   end
end
