classdef linearSolver < handle
   properties (Access = private)

      % Flag for Chronos existance
      ChronosFlag

      % Preconditioner Type
      PrecType

      % Preconditioner
      precOpt
      MfunL
      MfunR
      
      % Solver Type
      SolverType

      % Discretizer
      domain
      
      % Physics
      phys

      % Flag to request Preconditioner computation
      requestPrecComp
      Prec

      % starting vector
      x0

      % Statistics
      whenComputed
      aTimeComp
      aTimeSolve
      nSolve
      nComp
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
               fprintf('multiPhysics not yet supported\nFall back to matlab solver\n');
               obj.ChronosFlag = false;
               return
            end

            obj.domain = domainin;

            physname = obj.domain.solverNames(1);
            if(physname == 'SinglePhaseFlow' || physname == 'VariablySaturatedFlow' || physname == 'Poisson')
               obj.phys = 0;
            elseif(physname == 'Poromechanics')
               obj.phys = 1;
            else
               warning('Non supported Physics for linsolver, falling back to matlab solver');
               obj.ChronosFlag = false;
               return
            end

            % Chronos exists
            obj.ChronosFlag = true;
            addpath(genpath(ChronosDir));
            RACPDir = fullfile(gres_root,'..','aspamg_matlab','composed_precs');
            addpath(genpath(RACPDir));

            % Read XML
            if nargin > 1
               % Use input values
               data = readstruct(varargin{1},AttributeSuffix="");
            else
               % Get default values
               if obj.phys == 0
                  chronos_xml_default = fullfile(gres_root,'Code','linearSolver','chronos_xml_setup_CFD.xml');
               else
                  chronos_xml_default = fullfile(gres_root,'Code','linearSolver','chronos_xml_setup.xml');
               end

               % Read Defaults
               data = readstruct(chronos_xml_default,AttributeSuffix="");
            end
            
            % null or zero value to the unsued stuff
            obj.Prec   = [];
            obj.MfunL  = [];
            obj.MfunR  = [];
            obj.nSolve = 0;
            obj.nComp  = 0;
            obj.aTimeSolve    = 0;
            obj.aTimeComp     = 0;
            obj.whenComputed  = [];
            obj.precOpt       = -1;

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

            obj.params.tol   = data.general.tol;
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
            obj.requestPrecComp = true;
            obj.params.iter = -1;
            obj.params.lastRelres = obj.params.tol;
            obj.x0 = [];
         else
            obj.ChronosFlag = false;
         end
      end

      function printStats(obj)
         fprintf('Average Preconditioner computation time = %e\n',(obj.aTimeComp/obj.nComp));
         fprintf('Average Solve time = %e\n',(obj.aTimeSolve/obj.nSolve));
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

   end

end


