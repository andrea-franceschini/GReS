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
   end

   properties (Access = public)

      % Params struct
      params
   end

   methods (Access = public)
      
      % Constructor Function
      function obj = linearSolver(filename,domainin)
         
         % Check if chronos is available
         ChronosDir = fullfile(gres_root,'..','aspamg_matlab','sources');

         if isfolder(ChronosDir)
            
            if(numel(domainin.fields) > 1)
               fprintf('multiPhysics not yet supported\nFall back to matlab solver\n');
               obj.ChronosFlag = false;
               return
            end

            obj.domain = domainin;

            % Read XML
            data = readstruct(filename);
            
            % Chronos exists
            obj.ChronosFlag = true;
            addpath(genpath(ChronosDir));

            % null value to the preconditioner for now,
            obj.Prec = [];
            obj.MfunL = [];
            obj.MfunR = [];

            % Get the solver type
            obj.SolverType = lower(data.solver);
            
            % if GMRES get restart value
            if (obj.SolverType == "gmres")
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
                  obj.params.prolong    = data.prolong;
                  obj.params.coarsen   = data.coarsen;
                  obj.params.tspace   = data.tspace;
                  obj.params.filter   = data.filter;

               case 'fsai'
                  obj.params.smoother = data.smoother;
            end
         else
            obj.ChronosFlag = false;
         end
      end

      % Function for the computation of the preconditioner
      function computePrec(obj,A)

         % Chronos does not exist, use / to solve
         if ~obj.ChronosFlag
            return
         end

         if (norm(A-A',"fro")/norm(A,"fro") > 1e-7)
            obj.params.symm = false;
         else
            obj.params.symm = true;
         end

         switch obj.PrecType

            % Compute the AMG preconditioner
            case 'amg'
               
               % Treat Boundary conditions 
               lmax = eigs(A,1,'lm','FailureTreatment','keep','Display',0,'Tolerance',0.001,'MaxIterations',3);
               A(A==1) = lmax/10;

               % Compute the test space
               if(obj.domain.fields(1) == 'SinglePhaseFlow')
                  TV0 = ones(size(A,1),1);
               elseif(obj.domain.fields(1) == 'Poromechanics')
                  TV0 = mk_rbm_3d(obj.domain.grid.topology.coordinates);
               end
               
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

               % Actually compute the AMG
               obj.Prec = cpt_aspAMG(obj.params,A,TV0);

               % Define Mfun
               obj.MfunL = @(r) AMG_Vcycle(obj.Prec,A,r);
               obj.MfunR = @(r) r;

            % Compute the FSAI preconditioner
            case 'fsai'
               smootherOp = smoother(A,obj.params.smoother);

               % Define Mfun
               omega = smootherOp.omega;
               if strcmp(lower(obj.params.smoother.method),'afsai_enh')
                  F     = smootherOp.left;
                  FT    = smootherOp.right;
                  W     = smootherOp.W;
                  THETA = smootherOp.THETA;
                  obj.MfunL = @(x) omega*(FT*(F*x) + W*(THETA*(W'*x)));
                  obj.MfunR = @(x) omega*(FT*(F*x) + W*(THETA*(W'*x)));
               else
                  if smootherOp.LS_deg > 0
                     obj.MfunL = smootherOp.polyPrec;
                     obj.MfunR = smootherOp.polyPrec;
                  else
                     if strcmp(lower(obj.params.smoother.method),'blk_j')
                        obj.MfunL = smootherOp.BLKJ;
                        obj.MfunR = smootherOp.BLKJ;
                     elseif strcmp(lower(obj.params.smoother.method),'bafsai')
                        obj.MfunL = smootherOp.BAFSAI;
                        obj.MfunR = smootherOp.BAFSAI;
                     elseif strcmp(lower(obj.params.smoother.method),'ddsw')
                        obj.MfunL = smootherOp.DDSW1;
                        obj.MfunR = smootherOp.DDSW2;
                     else
                        if (numel(smootherOp.left_out) + numel(smootherOp.right_out)) == 0
                           % Simple smoother
                           obj.MfunL = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.left*x));
                           obj.MfunR = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.left*x));
                        else
                           % Nested smoother
                           obj.MfunL = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.right_out*...
                                               (smootherOp.left_out*(smootherOp.left*x))));
                           obj.MfunR = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.right_out*...
                                               (smootherOp.left_out*(smootherOp.left*x))));
                        end
                     end
                  end
               end
         end
      end


      % Function for the solution of the system
      % note that the A passed here might be slightly different than the one passed in the computation
      % if the two As differ too much the preconditioner loses effectiveness. Must be recomputed
      function [x,flag] = Solve(obj,A,b,nonLinIter)
         
         % Chronos does not exist, continue with matlab default
         if ~obj.ChronosFlag
            fprintf('hello from matlab solution\n');
            x = A\b;
            flag = 0;
            return
         end

         % Save the solver type
         firstSolver = obj.SolverType;

         % If not first iter the case is for sure nonlinear and for sure nonsymmetric
         if(nonLinIter > 1)
            obj.SolverType = 'gmres';
         else
            % If the matrix is nonSymmetric the use always GMRES
            if (norm(A-A',"fro")/norm(A,"fro") > 1e-7)
               obj.SolverType = 'gmres';
            end
         end

         switch obj.SolverType
            case 'gmres'

               % Solve the system by GMRES
               [x,flag,relres,obj.params.iter,resvec] = gmres(A,b,obj.params.restart,obj.params.tol,obj.params.maxit,obj.MfunL,obj.MfunR);

            case 'sqmr'

               % Solve the system by SQMR
               Afun = @(x) A*x;
               [x,flag,relres,obj.params.iter,resvec] = SQMR(Afun,b,obj.params.tol,obj.params.maxit,obj.MfunL,obj.MfunR);

         end

         % Reset the solver
         obj.SolverType = firstSolver;

      end
   end
end


