classdef linearSolver < handle
   properties (Access = private)

      % Flag for Chronos existance
      ChronosFlag

      % Preconditioner Type
      PrecType

      % Preconditioner
      Prec
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
            
            if(numel(domainin.fields) > 1)
               fprintf('multiPhysics not yet supported\nFall back to matlab solver\n');
               obj.ChronosFlag = false;
               return
            end

            obj.domain = domainin;

            if(obj.domain.fields(1) == 'SinglePhaseFlow' || obj.domain.fields(1) == 'VariablySaturatedFlow')
               obj.phys = 0;
            elseif(obj.domain.fields(1) == 'Poromechanics')
               obj.phys = 1;
            else
               warning('Non supported Physics for linsolver, falling back to matlab solver');
               obj.ChronosFlag = false;
               return
            end

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
            
            % Chronos exists
            obj.ChronosFlag = true;
            addpath(genpath(ChronosDir));

            % null value to the unsued stuff
            obj.Prec   = [];
            obj.MfunL  = [];
            obj.MfunR  = [];
            obj.nSolve = 0;
            obj.nComp  = 0;
            obj.aTimeSolve    = 0;
            obj.aTimeComp     = 0;
            obj.whenComputed  = [];

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

      % Function for the solution of the system
      % note that the A passed here might be slightly different than the one passed in the computation
      % if the two As differ too much the preconditioner loses effectiveness. Must be recomputed
      function [x,flag] = Solve(obj,A,b,time)
         
         % Chronos does not exist, continue with matlab default
         if ~obj.ChronosFlag || size(A,1) < 1e4 
            startT = tic;
            % Solve the system
            x = A\b;
            Tend = toc(startT);
            obj.aTimeSolve = obj.aTimeSolve + Tend;
            obj.nSolve = obj.nSolve + 1;
            flag = 0;
            obj.params.iter = 0;
            return
         end

         % Have the linear solver compute the Preconditioner if necessary
         if(obj.requestPrecComp || obj.params.iter > 500 || obj.params.lastRelres > obj.params.tol*1e3)
            obj.computePrec(A);
            obj.whenComputed(length(obj.whenComputed) + 1) = time;
            obj.params.iterSinceLastPrecComp = 0;
         else
            obj.params.iterSinceLastPrecComp = obj.params.iterSinceLastPrecComp + 1;
         end

         % Save the solver type
         firstSolver = obj.SolverType;

         % % If not first iter the case is for sure nonlinear and for sure nonsymmetric
         % if(nonLinIter > 1)
         %    obj.SolverType = 'gmres';
         % else
            % If the matrix is nonSymmetric the use always GMRES
            if (norm(A-A',"fro")/norm(A,"fro") > 1e-7)
               obj.SolverType = 'gmres';
            end
         % end

         startT = tic;
         switch obj.SolverType
            case 'gmres'

               % Solve the system by GMRES
               [x,flag,obj.params.lastRelres,iter1,resvec] = gmres_LEFT(A,b,obj.params.restart,obj.params.tol,...
                                                                        obj.params.maxit/obj.params.restart,obj.MfunL,obj.MfunR,obj.x0);
               obj.params.iter = (iter1(1) - 1) * obj.params.restart + iter1(2);

            case 'sqmr'

               % Solve the system by SQMR
               Afun = @(x) A*x;
               [x,flag,obj.params.lastRelres,obj.params.iter,resvec] = SQMR(Afun,b,obj.params.tol,obj.params.maxit,obj.MfunL,obj.MfunR,obj.x0);
         end

         Tend = toc(startT);
         obj.aTimeSolve = obj.aTimeSolve + Tend;
         obj.nSolve = obj.nSolve + 1;

         % Did not converge, if prec not computed for it try again
         if(flag == 1 && obj.params.iterSinceLastPrecComp > 0)
            fprintf('Trying to recompute the preconditioner to see if it manages to converge\n');
            obj.computePrec(A);
            obj.params.iterSinceLastPrecComp = 0;
            [x,flag] = obj.Solve(A,b,-1);
         end

         % Interesting problem
         if(flag == 1)
            x0 = obj.x0;
            fprintf('Iterations since last preconditioner computation %d\n',obj.params.iterSinceLastPrecComp);
            save('new_problem.mat','A','b','x0');
            error('Interesting problem spotted');
         end

         % Reset the solver
         obj.SolverType = firstSolver;

         % If the preconditioner has just been computed then do not compute it for the next iter
         if(obj.requestPrecComp)
            % Keep in memory the number of iter it did with the correct matrix
            obj.params.firstIterAfterPrecComp = obj.params.iter;
            obj.requestPrecComp = false;
         else
            % If the number of iterations changes too much then recompute the preconditioner
            if(obj.params.iter > 1.2 * obj.params.firstIterAfterPrecComp && obj.params.iter > obj.params.minIter)
               obj.requestPrecComp = true;
            end
         end

         % Store the new starting vector
         obj.x0 = x;
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
   end

   methods (Access = private)
      % Function for the computation of the preconditioner
      function computePrec(obj,A)

         if (norm(A-A',"fro")/norm(A,"fro") > 1e-7)
            obj.params.symm = false;
         else
            obj.params.symm = true;
         end

         time_start = tic;
         switch obj.PrecType

            % Compute the AMG preconditioner
            case 'amg'
               
               % Treat Boundary conditions 
               lmax = eigs(A,1,'lm','FailureTreatment','keep','Display',0,'Tolerance',0.001,'MaxIterations',3);
               A(A==1) = lmax/10;

               % Compute the test space
               if(obj.phys == 0) % fluids
                  TV0 = ones(size(A,1),1);
               elseif(obj.phys == 1)
                  TV0 = mk_rbm_3d(obj.domain.grid.topology.coordinates);
               end
               
               obj.set_DEBINFO();

               % Actually compute the AMG
               obj.Prec = cpt_aspAMG(obj.params,A,TV0);

               % Define Mfun
               obj.MfunL = @(r) AMG_Vcycle(obj.Prec,A,r);
               obj.MfunR = @(r) r;

            % Compute the FSAI preconditioner
            case 'fsai'
               smootherOp = smoother(A,obj.params.symm,obj.params.smoother);

               % Define Mfun
               [obj.MfunL,obj.MfunR] = obj.defineMfunFSAI(smootherOp);
         end
         T_setup = toc(time_start);
         obj.aTimeComp = obj.aTimeComp + T_setup;
         obj.nComp = obj.nComp + 1;
      end

      % Function to determine how MfunL and MfunR are for each fsai preconditioner
      function [MfunL,MfunR] = defineMfunFSAI(obj,smootherOp)
         omega = smootherOp.omega;
         if strcmp(lower(obj.params.smoother.method),'afsai_enh')
            F     = smootherOp.left;
            FT    = smootherOp.right;
            W     = smootherOp.W;
            THETA = smootherOp.THETA;
            MfunL = @(x) omega*(FT*(F*x) + W*(THETA*(W'*x)));
            MfunR = @(x) omega*(FT*(F*x) + W*(THETA*(W'*x)));
         else
            if smootherOp.LS_deg > 0
               MfunL = smootherOp.polyPrec;
               MfunR = smootherOp.polyPrec;
            else
               if strcmp(lower(obj.params.smoother.method),'blk_j')
                  MfunL = smootherOp.BLKJ;
                  MfunR = smootherOp.BLKJ;
               elseif strcmp(lower(obj.params.smoother.method),'bafsai')
                  MfunL = smootherOp.BAFSAI;
                  MfunR = smootherOp.BAFSAI;
               elseif strcmp(lower(obj.params.smoother.method),'ddsw')
                  MfunL = smootherOp.DDSW1;
                  MfunR = smootherOp.DDSW2;
               else
                  if (numel(smootherOp.left_out) + numel(smootherOp.right_out)) == 0
                     % Simple smoother
                     MfunL = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.left*x));
                     MfunR = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.left*x));
                  else
                     % Nested smoother
                     MfunL = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.right_out*...
                                         (smootherOp.left_out*(smootherOp.left*x))));
                     MfunR = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.right_out*...
                                         (smootherOp.left_out*(smootherOp.left*x))));
                  end
               end
            end
         end
      end
   end

   methods (Static)
      function set_DEBINFO(obj)
         global DEBINFO;
         % PARTE GENERALE
         DEBINFO.flag = true;
         DEBINFO.flag = false;
         % PARTE PER LA PROLONGATION
         DEBINFO.prol = [];
         % STAMPARE SI/NO
         DEBINFO.prol.prt_flag = false;
         % UNITA DI STAMPA
         DEBINFO.prol.ofile = 0;
         % STAMPA NUMERO DI ITERAZIONI NEL CALCOLO PROL
         DEBINFO.prol.it_print = false;
         % STAMPA LISTA VICINI NEL CALCOLO PROL
         DEBINFO.prol.neigh_print = false;
         
         % PARTE PER IL COARSENING
         DEBINFO.coarsen = [];
         % STAMPARE SI/NO
         DEBINFO.coarsen.draw_dist = false;
      end
   end
end


