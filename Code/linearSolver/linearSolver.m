classdef linearSolver < handle
   properties (Access = private)

      % Flag for Chronos existance
      ChronosFlag

      % Preconditioner Type
      PrecType

      % Preconditioner
      Prec
      
      % Solver Type
      SolverType
   end

   properties (Access = public)

      % Params struct
      params
   end

   methods (Access = public)
      
      % Constructor Function
      function obj = linearSolver(filename)
         
         % Check if chronos is available
         ChronosDir = fullfile(gres_root,'..','aspAMGciao_matlab');

         if isfolder(targetFolder)
            
            data = readstruct(filename);
            
            % Chronos exists
            ChronosFlag = true;
            addpath(genpath(ChronosDir));

            % null value to the preconditioner for now,
            obj.Prec = [];

            % Get the solver type
            SolverType = lower(data.solver);
            
            % if GMRES get restart value
            if (SolverType == 'gmres')
               params.restart = params.general.restart;
            end

            % Get the preconditioner type
            PrecType = lower(data.preconditioner);

            params.tol   = data.general.tol;
            params.maxit = data.general.maxit;

            % Get the different parameters according to the prectype
            switch PrecType
               case 'amg'
                  params.amg      = data.amg;
                  params.smoother = data.smoother;
                  params.prolo    = data.prolong;
                  params.coarse   = data.coarse;
                  params.tspace   = data.tspace;
                  params.filter   = data.filter;

               case 'fsai'
                  params.smoother = data.smoother;
               end
            end
         else
            ChronosFlag = false;
         end
      end

      % Function for the computation of the preconditioner
      function computePrec(A,varargin)

         % Chronos does not exist, use / to solve
         if ~ChronosFlag
            return
         end

         switch PrecType
            % Compute the AMG preconditioner
            case 'amg'
               
               if varargin < 1
                  error('linearSolver::computePrec','necessary coordinates to compute test space');
               end

               % Treat Boundary conditions 
               lmax = eigs(A,1,'lm','FailureTreatment','keep','Display',0,'Tolerance',0.001,'MaxIterations',3);
               A(A==1) = lmax/10;

               % Compute the test space
               TV0 = mk_rbm_3d(varargin{1});
               
               % Actually compute the AMG
               Prec = cpt_aspAMG(param,A,TV0);

            % Compute the FSAI preconditioner
            case 'fsai'
               Prec = smoother(A,param);
            end 
         end
      end


      % Function for the solution of the system
      % note that the A passed here might be slightly different than the one passed in the computation
      % if the two As differ too much the preconditioner loses effectiveness. Must be recomputed
      function [flag] = Solve(A,x,b,nonLinIter)
         
         % Chronos does not exist, continue with matlab default
         if ~ChronosFlag
            x = A\b;
            flag = 0;
            return
         end

         % Save the solver type
         firstSolver = SolverType;

         % If not first iter the case is for sure nonlinear and for sure nonsymmetric
         if(nonLinIter > 1)
            SolverType = 'gmres';
         else
            % If the matrix is nonSymmetric the use always GMRES
            if (norm(A-A',"fro")/norm(A,"fro") > 1e-7)
               SolverType = 'gmres';
            end
         end

         switch SolverType
            case 'gmres'

               % Solve the system by GMRES
               Mfun = @(r) AMG_Vcycle(AMG_prec,A,r);
               [x,flag,relres,param.iter,resvec] = gmres(A,b,params.restart,params.tol,params.itmax,Mfun);

            case 'sqmr'

               % Solve the system by SQMR
               Afun = @(x) A*x;
               Mfun = @(r) AMG_Vcycle(AMG_prec,A,r);
               IDfun = @(x) x;
               [x,flag,relres,param.iter,resvec] = SQMR(Afun,b,params.tol,params.itmax,Mfun,IDfun);

         end

         % Reset the solver
         SolverType = firstSolver;

      end
   end
end

