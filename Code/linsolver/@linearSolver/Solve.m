
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
   if ~obj.ChronosFlag || (getGlobalSize(A) < obj.matlabMaxSize)
      [x,flag] = matlab_solve(obj,A,b);
      return
   end

   % Contact has opened a fracture or something similar so amg does not converge. 
   % Directly recompute the preconditioner
   if obj.Prec.phys == 1.1 
      if obj.generalsolver.iterConfig > obj.iterConfigOld
         obj.requestPrecComp = true;
         obj.iterConfigOld = obj.generalsolver.iterConfig;
      elseif obj.generalsolver.iterConfig < obj.iterConfigOld
         obj.iterConfigOld = obj.generalsolver.iterConfig;
      end
   end

   % Save the solver type
   firstSolver = obj.SolverType;

   [globalsymm,maxval,symMat] = checkSymmetry(A,obj.nsyTol);
   
   if globalsymm == 0
      % If the matrix is nonSymmetric the use always GMRES
      obj.SolverType = 'gmres';
      if obj.DEBUGflag
         fprintf("The matrix is nonsymmetric with a maximum nonsymmetry of %e\n",maxval);
      end
   end

   % Have the linear solver compute the Preconditioner if necessary
   if(obj.requestPrecComp || obj.params.iter > 600 || obj.params.lastRelres > obj.params.tol*1e3)
      if obj.DEBUGflag
         fprintf("Computing the preconditioner\n");
      end

      time_start = tic;
      obj.Prec.Compute(A,symMat);
      T_setup = toc(time_start);

      obj.aTimeComp = obj.aTimeComp + T_setup;
      obj.nComp = obj.nComp + 1;
      obj.whenComputed(length(obj.whenComputed) + 1) = time;
      obj.params.iterSinceLastPrecComp = 0;
   else
      obj.params.iterSinceLastPrecComp = obj.params.iterSinceLastPrecComp + 1;
   end

   %save('Ab.mat',"A","b");
   if iscell(A)
      Amat = cell2matrix(A);
   end

   startT = tic;
   switch obj.SolverType
      case 'gmres'

         % Solve the system by GMRES
         [x,flag,obj.params.lastRelres,iter1,resvec] = gmres_RIGHT(Amat,b,obj.params.restart,obj.params.tol,...
                                                                   obj.params.maxit/obj.params.restart,...
                                                                   obj.Prec.Apply_L,obj.Prec.Apply_R,obj.x0,obj.DEBUGflag);
         obj.params.iter = (iter1(1) - 1) * obj.params.restart + iter1(2);
         
      case 'sqmr'

         % Solve the system by SQMR
         Afun = @(x) Amat*x;
         [x,flag,obj.params.lastRelres,obj.params.iter,resvec] = SQMR(Afun,b,obj.params.tol,obj.params.maxit,...
                                                                      obj.Prec.Apply_L,obj.Prec.Apply_R,obj.x0,obj.DEBUGflag);

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
      return;
   end

   % Interesting problem
   if(flag == 1)
      if obj.DEBUGflag
         fprintf('Iterations since last preconditioner computation %d\n',obj.params.iterSinceLastPrecComp);
      end
      [~,~] = matlab_solve(obj,A,b);
      TV0 = obj.Prec.TV0;
      save('new_problem.mat','A','b','TV0');
      error('matlab could solve and chronos did not.');
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
      if(obj.params.iter > 1.5 * obj.params.firstIterAfterPrecComp && obj.params.iter > obj.params.minIter)
         obj.requestPrecComp = true;
      end
   end

   % Store the new starting vector
   obj.x0 = x;
end




















% Function to get the global size of the system
function [gSize] = getGlobalSize(A)

   [nBlockRows,nBlockCols] = size(A);
   rowSizes = zeros(nBlockRows,1);
   
   for i = 1:nBlockRows
      for j = 1:nBlockCols
         C = A{i,j};
         if ~isempty(C)
            rowSizes(i) = size(C,1);
            break;
         end
      end
   end

   gSize = sum(rowSizes);
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

   if obj.DEBUGflag
      fprintf('condition number of the matrix %e\n',condest(A));
   end

   obj.aTimeSolve = obj.aTimeSolve + Tend;
   obj.nSolve = obj.nSolve + 1;
   flag = 0;
   obj.params.iter = 0;
end

function [globalsymm,maxval,symMat] = checkSymmetry(A,eps1)
   
   % Base Case: Numeric Matrix
   if ~iscell(A)
       diff_mat = abs(A - A') - eps1 .* max(abs(A), abs(A'));
       diff_vec = diff_mat(diff_mat > 0);
       
       if isempty(diff_vec)
           maxval = 0;
           globalsymm = 1;
       else
           maxval = max(diff_vec);
           globalsymm = 0;
       end
       symMat = globalsymm;
       return
   end
   
   % Allocate the stuff
   N = size(A,1);
   symm = zeros(sum(1:N),1);
   val = zeros(sum(1:N),1);
   cont = 1;
   
   % Loop over the blocks
   for j = 1:N
      for i = 1:j
         if i == j
            % Diagonal Block
            [symm(cont), val(cont)] = checkSymmetry(A{i,i},eps1);
         else
            % Off-Diagonal Block
            diff_mat = abs(A{i,j} - A{j,i}') - eps1 .* max(abs(A{i,j}), abs(A{j,i}'));
            diff_vec = diff_mat(diff_mat > 0);
            
            if isempty(diff_vec)
                symm(cont) = 1;
                val(cont) = 0;
            else
                symm(cont) = 0;
                val(cont) = max(diff_vec);
            end
         end
         cont = cont + 1;
      end
   end
   
   % Global mapping and output
   symMat = zeros(N, N);
   symMat(triu(true(N))) = symm;
   symMat = symMat + triu(symMat, 1).';
   
   globalsymm = min(symm);
   maxval = max(val);
end