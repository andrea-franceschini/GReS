classdef linearSolver < handle
   properties (Access = private)

      % Flag for debug
      DEBUGflag = true
      matlabMaxSize = 1e1
      nsyTol = 1e-15

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

      % Params struct
      params
   end

   methods (Access = public)

      % Constructor Function
      function obj = linearSolver(generalsolver,usrInput,physname)

        [isAMGavailable,AMGdir] = obj.isAMGavailable();

         if isAMGavailable

            obj.generalsolver = generalsolver;

            % Create the preconditioner object, check if the physics is supported
            [obj.Prec,obj.ChronosFlag] = preconditioner.create(obj.DEBUGflag,obj.nsyTol,generalsolver,usrInput,physname);

            % Non supported physics for the preconditioner
            if ~obj.ChronosFlag
               return;
            end

            % Chronos exists
            addpath(genpath(AMGdir));
            RACPDir = fullfile(gres_root,'ThirdPartyLibs','aspAMG','composed_precs','RACP');
            addpath(genpath(RACPDir));

            % First time solving request preconditioner computation
            obj.params.iter = -1;
            obj.params.lastRelres = 1e10;
            obj.params.tol   = generalsolver.simparams.relTol;

            % Get default values
            chronos_xml_default = fullfile(gres_root,'Code','linsolver','XML_setup','chronos_xml_setup.xml');

            % Read Defaults
            data = readstruct(chronos_xml_default,AttributeSuffix="");

            % Get the solver type
            obj.SolverType = lower(data.solver);
            obj.params.maxit = data.general.maxit;
            obj.params.minIter = obj.Prec.params.minIter;

            % if GMRES get restart value
            if (obj.SolverType == "gmres")
               obj.params.restart = data.general.restart;
            else
               obj.params.restart = 100;
            end
         end
      end

      function printStats(obj)
         fprintf('Used %d threads during mex\n',obj.Prec.maxThreads);
         fprintf('Average Preconditioner computation time = %e\n',(obj.aTimeComp/obj.nComp));
         fprintf('Average Solve time = %e\n',(obj.aTimeSolve/obj.nSolve));
         fprintf('Average number of iterations = %f\n',(obj.aIter/obj.nSolve));
         fprintf('Max number of iterations = %d\n',obj.maxIter);
         fprintf('The preconditioner was computed at time(s):\n');
         for i = 1:length(obj.whenComputed)
            fprintf('             %d\t%e\n',i,obj.whenComputed(i));
         end
         fprintf('Total time for computation of the linear systems = %e\n',obj.aTimeComp+obj.aTimeSolve);
      end

      % Function to solve the system
      [x,flag] = Solve(obj,A,b,time)

   end


   methods (Static)
     function [out,varargout] = isAMGavailable()
       % Check if aspAMG is available in the TPL
       aspAMGdir = fullfile(gres_root,'ThirdPartyLibs','aspAMG','sources');
       out = isfolder(aspAMGdir);

       if nargout > 1
         varargout{1} = aspAMGdir;
       end
     end
   end
end


