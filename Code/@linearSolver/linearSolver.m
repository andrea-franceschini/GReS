classdef linearSolver < handle
   properties (Access = private)

      % Flag for debug
      DEBUGflag = true
      matlabMaxSize = 5e4
      nsyTol = 1e-12

      % Flag for Chronos existance
      ChronosFlag = false

      % Preconditioner Type
      PrecType

      % Preconditioner
      Prec = []

      % Preconditioner application
      MfunL = []
      MfunR = []
      
      % Flag to request Preconditioner computation
      requestPrecComp = true

      % starting vector
      x0 = []

      % vector containing the insufficient tolerances
      notSuffTol = []

      % Solver Type
      SolverType

      % Discretizer
      domain

      % Number of domains/interfaces
      nDom
      nInt
      
      % Physics
      phys
      
      % Flag to treat multiple domains as multiple domains
      multidomFlag = false

      % Flag to know if the problem has multiphysics
      multiPhysFlag = false

      % Statistics
      whenComputed = []
      aTimeComp = 0
      aTimeSolve = 0
      nSolve = 0
      nComp = 0
      maxIter = -1
      aIter = 0

      % Max Threads
      maxThreads
   end


   properties (Access = public)

      % Params struct
      params
   end

   methods (Access = public)
      
      % Constructor Function
      function obj = linearSolver(domainin,varargin)
         
         % Check if chronos is available
         ChronosDir = fullfile(gres_root,'..','aspamg_matlab','sources');

         if isfolder(ChronosDir)
            
            if(numel(domainin.physicsSolvers) > 1)

               obj.multiPhysFlag = true;

               if obj.DEBUGflag
                  fprintf('multiPhysics not yet supported\nFall back to matlab solver\n');
               end
               return
            end

            if nargin >= 2
               interfacein = varargin{1};
            else
               interfacein = ones(1);
            end

            obj.domain = domainin;
            obj.nDom = length(domainin);
            obj.nInt = length(interfacein);

            physname = obj.domain(1).solverNames(1);
            if(contains(physname,'SinglePhaseFlow') || physname == 'VariablySaturatedFlow' || physname == 'Poisson')
               obj.phys = 0;
            elseif(physname == 'Poromechanics')
               obj.phys = 1;
            else
               if obj.DEBUGflag
                  warning('Non supported Physics for linsolver, falling back to matlab solver');
                  fprintf('physics: %s\n',physname);
               end
               return
            end

            % Chronos exists
            obj.ChronosFlag = true;
            addpath(genpath(ChronosDir));
            RACPDir = fullfile(gres_root,'..','aspamg_matlab','composed_precs','RACP');
            addpath(genpath(RACPDir));

            % Read XML
            if nargin > 2
               % Use input values
               data = readstruct(varargin{2},AttributeSuffix="");
            else
               % Get default values
               if obj.phys == 0
                  chronos_xml_default = fullfile(gres_root,'Code','@linearSolver','XML_setup','chronos_xml_setup_CFD.xml');
               else
                  chronos_xml_default = fullfile(gres_root,'Code','@linearSolver','XML_setup','chronos_xml_setup.xml');
               end

               % Read Defaults
               data = readstruct(chronos_xml_default,AttributeSuffix="");
            end
            
            % Get the solver type
            obj.SolverType = lower(data.solver);
            
            % if GMRES get restart value
            if (obj.SolverType == 'gmres')
               obj.params.restart = data.general.restart;
            else
               obj.params.restart = 100;
            end

            % Get the preconditioner type
            obj.PrecType = lower(data.preconditioner);

            obj.params.tol   = obj.domain.simparams.relTol;
            obj.params.maxit = data.general.maxit;

            % Get the different parameters according to the prectype
            switch obj.PrecType
               case 'amg'
                  obj.params.amg      = data.amg;
                  obj.params.smoother = data.smoother;
                  obj.params.prolong  = data.prolong;
                  obj.params.coarsen  = data.coarsen;
                  obj.params.tspace   = data.tspace;
                  obj.params.filter   = data.filter;
                  obj.params.minIter  = 30;

               case 'fsai'
                  obj.params.smoother = data.smoother;
                  obj.params.minIter  = 300;
            end

            % First time solving request preconditioner computation
            obj.params.iter = -1;
            obj.params.lastRelres = 1e10;
   
            % Set maximum number of threads to use if the system provides less
            obj.maxThreads = maxNumCompThreads;
            obj.params.smoother.nthread = min(obj.params.smoother.nthread,obj.maxThreads);
            obj.params.prolong.np = min(obj.params.prolong.np,obj.maxThreads);
            obj.params.filter.np = min(obj.params.filter.np,obj.maxThreads);

         end
      end

      function printStats(obj)
         fprintf('Average Preconditioner computation time = %e\n',(obj.aTimeComp/obj.nComp));
         fprintf('Average Solve time = %e\n',(obj.aTimeSolve/obj.nSolve));
         fprintf('Average number of iterations = %e\n',(obj.aIter/obj.nSolve));
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

   methods (Access = private)

      % Function to compute the preconditioner
      computePrec(obj,A)

      % Function to compute the preconditioner for the single block (single physics)
      computeSinglePhPrec(obj,A);

      % Function to compute the RACP preconditioner for the lagrange multiplier case (single physics multi domain)
      computeRACP(obj,A)
   end

end


