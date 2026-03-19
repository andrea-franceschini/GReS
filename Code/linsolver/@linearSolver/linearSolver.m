classdef linearSolver < handle
   properties (SetAccess = private, GetAccess = public)

      % Flag for debug
      DEBUGflag = false
      matlabMaxSize = 1e5
      nsyTol = 100*eps
      fullInfo = false

      % Flag for Chronos existance
      ChronosFlag = false

      % Flag to request Preconditioner computation
      requestPrecComp = true

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
      
      % Full info stats
      iterLin = []
      solveTLin = []
      symFlagLin = []
      precCompLin = []
      newtonLin = []
      timeLin = []

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
               fprintf('The user requested the use of matlab\n');
               return;
            end
         end

         if isfolder(ChronosDir)
            if strcmp(computer('arch'),'maca64')
               fileMex = fullfile(ChronosDir,'Preconditioner','AMG','filter','MEX_Prol_Filter','FilterProl_wrap.mexmaca64');
            else
               fileMex = fullfile(ChronosDir,'Preconditioner','AMG','filter','MEX_Prol_Filter','FilterProl_wrap.mexa64');
            end
            
            if ~isfile(fileMex)
               warning('Chronos_Lab submodule is present, but not compiled. Using matlab fallback');
               return;
            end

            obj.generalsolver = generalsolver;

            % Create the preconditioner object, check if the physics is supported
            [obj.Prec,obj.ChronosFlag] = preconditioner.create(obj.DEBUGflag,obj.nsyTol,generalsolver,physname);

            % Non supported physics for the preconditioner
            if ~obj.ChronosFlag
               return;
            end

            % Chronos exists
            addpath(genpath(ChronosDir));

            % First time solving request preconditioner computation
            obj.params.iter = -1;
            obj.params.lastRelres = 1e10;
            obj.params.tol = generalsolver.simparams.relTol;

            % Get default values
            chronos_xml_default = fullfile(gres_root,'Code','linsolver','XML_setup','chronos_xml_setup.xml');

            % Read Defaults
            data = readInput(chronos_xml_default);

            % Get the solver type
            obj.SolverType = lower(data.solver);
            obj.params.maxit = data.general.maxit;
            obj.params.minIter = obj.Prec.params.minIter;

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

      function printStats(obj)
         fprintf('\n\n\nAverage Preconditioner computation time = %e\n',(obj.aTimeComp/obj.nComp));
         fprintf('The preconditioner was computed at time(s):\n');
         for i = 1:length(obj.whenComputed)
            fprintf('             %d\t%e\n',i,obj.whenComputed(i));
         end
         fprintf('\nAverage Solve time = %e\n',(obj.aTimeSolve/obj.nSolve));
         fprintf('\nAverage number of iterations = %f\n',(obj.aIter/obj.nSolve));
         fprintf('Max number of iterations = %d\n',obj.maxIter);
         fprintf('\nTotal time for computation of the linear systems = %e\n',obj.aTimeComp+obj.aTimeSolve);

         if obj.fullInfo
            fprintf('\nUsed %d threads during mex\n',obj.Prec.maxThreads);
            fprintf('\n-------------------------------------------------------------------------\n')
            fprintf('| %11s | %10s | %4s | %13s | %4s | %13s |\n','Time','NewtonIter','Iter','Time','Symm','PrecTime');
            fprintf('-------------------------------------------------------------------------\n')
            timeOld = obj.timeLin(1);
            for i = 1:size(obj.solveTLin,2)
               if timeOld ~= obj.timeLin(i)
                  fprintf('-------------------------------------------------------------------------\n')
               end
               fprintf('| %.5e | %10d | %4d | %.7e | %4d | %.7e |\n',obj.timeLin(i),obj.newtonLin(i),obj.iterLin(i),obj.solveTLin(i),obj.symFlagLin(i),obj.precCompLin(i));
            end
         end
      end

      % Function to solve the system
      [x,flag] = SolveLin(obj,A,b,time)

   end
end


