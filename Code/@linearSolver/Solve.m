
% Function for the solution of the system
% note that the A passed here might be slightly different than the one passed in the computationi of the preconditioner
% if the two As differ too much the preconditioner loses effectiveness. Must be recomputed

% A is passed directly as a cell array, meaning it is already split in the various blocks (A11,A12,A21,A22 for a 
% single physics single domain with lagrange multipliers) 
function [x,flag] = Solve(obj,A,b,time)
   
   if obj.DEBUGflag
      A
   end

   % Chronos does not exist, continue with matlab default
   if ~obj.ChronosFlag || size(A{1,1},1) < obj.matlabMaxSize
      [x,flag] = matlab_solve(obj,A,b);
      return
   end

   % Have the linear solver compute the Preconditioner if necessary
   if(obj.requestPrecComp || obj.params.iter > 600 || obj.params.lastRelres > obj.params.tol*1e3)
      obj.computePrec(A);
      obj.whenComputed(length(obj.whenComputed) + 1) = time;
      obj.params.iterSinceLastPrecComp = 0;
   else
      obj.params.iterSinceLastPrecComp = obj.params.iterSinceLastPrecComp + 1;
   end

   % Save the solver type
   firstSolver = obj.SolverType;

   %save('Ab.mat',"A","b");
   if iscell(A)
      Amat = cell2matrix(A);
   end

   % If the matrix is nonSymmetric the use always GMRES
   if (norm(Amat-Amat','f')/norm(Amat,'f') > obj.nsyTol)
      if obj.DEBUGflag
         fprintf('\nsym = %e\n\n',norm(Amat-Amat','f')/norm(Amat,'f'));
      end
      obj.SolverType = 'gmres';
   end

   startT = tic;
   switch obj.SolverType
      case 'gmres'

         % Solve the system by GMRES
         [x,flag,obj.params.lastRelres,iter1,resvec] = gmres_RIGHT(Amat,b,obj.params.restart,obj.params.tol,...
                                                                   obj.params.maxit/obj.params.restart,obj.MfunL,obj.MfunR,obj.x0,obj.DEBUGflag);
         obj.params.iter = (iter1(1) - 1) * obj.params.restart + iter1(2);
         
%          if obj.params.lastRelres > obj.params.tol
%             obj.notSuffTol(length(obj.notSuffTol)+1) = obj.params.tol/obj.params.lastRelres;
% 
%             [x,flag,obj.params.lastRelres,iter1,resvec] = gmres_RIGHT(Amat,b,obj.params.restart,obj.params.tol*obj.notSuffTol(end)*0.1,...
%                                                                       obj.params.maxit/obj.params.restart,obj.MfunL,obj.MfunR,x,obj.DEBUGflag);
%             obj.params.iter = obj.params.iter + (iter1(1) - 1) * obj.params.restart + iter1(2);
%          end


      case 'sqmr'

         % Solve the system by SQMR
         Afun = @(x) Amat*x;
         [x,flag,obj.params.lastRelres,obj.params.iter,resvec] = SQMR(Afun,b,obj.params.tol,obj.params.maxit,obj.MfunL,obj.MfunR,obj.x0,obj.DEBUGflag);

%          if obj.params.lastRelres > obj.params.tol
%             obj.notSuffTol(length(obj.notSuffTol)+1) = obj.params.tol/obj.params.lastRelres;
%             if obj.DEBUGflag
%                fprintf('recomputing solution\n');
%             end
%             
%             [x,flag,obj.params.lastRelres,obj.params.iter,resvec] = SQMR(Afun,b,obj.params.tol,obj.params.maxit,obj.MfunL,obj.MfunR,obj.x0,obj.DEBUGflag);
%          end
   end

   Tend = toc(startT);
   obj.aTimeSolve = obj.aTimeSolve + Tend;
   obj.nSolve = obj.nSolve + 1;
   obj.aIter = obj.aIter + obj.params.iter;
   obj.maxIter = max(obj.maxIter,obj.params.iter);

   % Did not converge, if prec not computed for it try again
   if(flag == 1 && obj.params.iterSinceLastPrecComp > 0)
      if obj.DEBUGflag
         fprintf('Trying to recompute the preconditioner to see if it manages to converge\n');
      end
      obj.params.iterSinceLastPrecComp = 0;
      obj.requestPrecComp = true;
      [x,flag] = obj.Solve(A,b,time);
   end

   % Interesting problem
   if(flag == 1)
      x0 = obj.x0;
      if obj.DEBUGflag
         fprintf('Iterations since last preconditioner computation %d\n',obj.params.iterSinceLastPrecComp);
      end
      %save('new_problem.mat','A','b','x0');
      error('');
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




function [x,flag] = matlab_solve(obj,A,b)

   if obj.DEBUGflag
      fprintf('Fallback to matlab due to size or chronos inexistance\n');
   end
   startT = tic;
   % Solve the system
   A = cell2matrix(A);
   x = A\b;
   Tend = toc(startT);
   obj.aTimeSolve = obj.aTimeSolve + Tend;
   obj.nSolve = obj.nSolve + 1;
   flag = 0;
   obj.params.iter = 0;
end
