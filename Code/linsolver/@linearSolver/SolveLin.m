function [x,flag] = SolveLin(obj,A,b,time)
% This file implements a Solve method and related utilities intended to be
% used as part of a linear solver class (linsolver). The Solve function
% orchestrates choosing between internal MATLAB direct solve and an external
% iterative solver (Chronos), manages preconditioner computation and reuse,
% gathers timing/iteration statistics, and handles fallback/retry logic.
%
% Relationship with linsolver class:
% - The code expects to be a method of an object "obj" that encapsulates the
%   solver state and configuration.
% - Prec is expected to be an object with methods Compute and function
%   handles Apply_L and Apply_R used as left/right preconditioners.
% - The linsolver class should provide gmres_RIGHT and SQMR wrappers or
%   have them available on the path. The code also relies on helper
%   utilities cell2matrix, fixPattern, and checkSymmetry.
% - Solve manages recursive calls to itself when it requests a
%   recomputation of the preconditioner and retrying the solve.
%
% Usage:
% - Call as: [x,flag] = obj.Solve(Acell,b,currentTime)
%   where Acell is a cell array of blocks representing the system matrix,
%   b is the right-hand side, and currentTime is a scalar timestamp used
%   for profiling (e.g., simulation time).
%
% Notes and behaviour details:
% - For small problems or when Chronos (external iterative solver) is not
%   available, the method falls back to matlab_solve which uses the direct
%   backslash on the assembled matrix.
% - If the matrix is detected non-symmetric, SolverType is forced to 'gmres'.
% - The preconditioner is computed when requested or when convergence
%   metrics indicate poor performance (high iteration counts or large
%   relative residual). The decision logic is encoded using fields in obj.
% - After a failed solve (flag == 1), the method may attempt to recompute
%   the preconditioner and retry once. If MATLAB direct solve succeeds but
%   Chronos does not, the code saves a snapshot ('new_problem.mat') and
%   raises an error to aid debugging.
% - Note that the A passed here might be slightly different than the one passed in the computation of the preconditioner
%   if the two As differ too much the preconditioner loses effectiveness. Must be recomputed
% - A is passed directly as a cell array, meaning it is already split in the various blocks (A11,A12,A21,A22 for a
%   single physics single domain with lagrange multipliers)
   
   if obj.DEBUGflag
      A
   end

   oldProbSize = size(obj.x0,1);
   sizeDiff = 0;

   % Chronos does not exist, continue with matlab default
   if ~obj.ChronosFlag || (getGlobalSize(A) < obj.matlabMaxSize)
      [x,flag] = matlab_solve(obj,A,b);
      return
   elseif oldProbSize ~= 0
      newProbSize = size(b,1);
      sizeDiff = newProbSize - oldProbSize;
      if sizeDiff > 0
         obj.x0 = [obj.x0; zeros(sizeDiff,1)];
      fprintf('changed size from %d to %d\n',oldProbSize,newProbSize);
      end
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

   % Apply Ruix Scaling on the Block matrix
   if obj.nIterRuiz > 0
      % Compute and apply Ritz scaling on A
      [A,obj.Prec.D] = ruiz_block_symmetric(A,obj.nIterRuiz,obj.tolRuiz,obj.DEBUGflag);

      % Prepare the D for future applications
      obj.Prec.D = cellfun(@(Di) diag(Di), obj.Prec.D, 'UniformOutput', false);

      % Single vector D to apply to the rhs
      D = diag(vertcat(obj.Prec.D{:}));

      % Apply the scaling to the rhs
      b = D*b;
   else
      obj.Prec.D = {};
   end

   % Fix the pattern to be symmetric and check the symmetry of the
   % resulting matrix
   [A] = fixPattern(A);
   [globalsymm,maxval,symMat] = checkSymmetry(A,obj.nsyTol);
   
   if globalsymm == 0
      % If the matrix is nonSymmetric then use always GMRES
      obj.SolverType = 'gmres';
      gresLog().log(3,'The matrix is nonsymmetric with a maximum nonsymmetry of %e\n',maxval);
   end

   % Have the linear solver compute the Preconditioner if necessary
   if(obj.requestPrecComp || obj.params.iter > 600 || obj.params.lastRelres > obj.params.tol*1e3)
      gresLog().log(3,'Computing the preconditioner\n');

      time_start = tic;
      obj.Prec.Compute(A,symMat);
      T_setup = toc(time_start);

      obj.aTimeComp = obj.aTimeComp + T_setup;
      obj.nComp = obj.nComp + 1;
      obj.whenComputed(length(obj.whenComputed) + 1) = time;
      obj.params.iterSinceLastPrecComp = 0;
      obj.sizeDiff = 0;
      gresLog().log(3,'Finished computing the preconditioner\n');
   else
      obj.params.iterSinceLastPrecComp = obj.params.iterSinceLastPrecComp + 1;
   end

   if iscell(A)
      Amat = cell2matrix(A);
      if obj.DEBUGflag
         symValue = norm(Amat-Amat','f')/norm(Amat,'f');
      end
   end

   if sizeDiff > 0 
      obj.sizeDiff = obj.sizeDiff + sizeDiff;
      
      obj.precL = @(x) [obj.Prec.Apply_L(x(1:end-obj.sizeDiff)); x(end-obj.sizeDiff+1:end)];
   elseif isempty(obj.precL)
      obj.precL = obj.Prec.Apply_L;
   end

   startT = tic;
   switch obj.SolverType
      case 'gmres'

         % Solve the system by GMRES
         [x,flag,obj.params.lastRelres,iter1,resvec] = gmres_RIGHT(Amat,b,obj.params.restart,obj.params.tol,...
                                                                   obj.params.maxit/obj.params.restart,...
                                                                   obj.precL,obj.Prec.Apply_R,obj.x0,obj.DEBUGflag);
         obj.params.iter = (iter1(1) - 1) * obj.params.restart + iter1(2);
         
      case 'sqmr'

         % Solve the system by SQMR
         Afun = @(x) Amat*x;
         [x,flag,obj.params.lastRelres,obj.params.iter,resvec] = SQMR(Afun,b,obj.params.tol,obj.params.maxit,...
                                                                      obj.precL,obj.Prec.Apply_R,obj.x0,obj.DEBUGflag);

   end

   % De apply ruiz from the result
   if obj.nIterRuiz > 0
      x = D*x;
   end

   Tend = toc(startT);

   % Save statistics for profiling or info in general
   obj.aTimeSolve = obj.aTimeSolve + Tend;
   obj.nSolve = obj.nSolve + 1;
   obj.aIter = obj.aIter + obj.params.iter;
   obj.maxIter = max(obj.maxIter,obj.params.iter);

   if obj.fullInfo
      obj.iterLin(obj.nSolve) = obj.params.iter;
      obj.timeLin(obj.nSolve) = time; 
      obj.solveTLin(obj.nSolve) = Tend;
      if exist('obj.generalsolver.iterNL','var')
         obj.newtonLin(obj.nSolve) = obj.generalsolver.iterNL;
      else
         obj.newtonLin(obj.nSolve) = obj.nSolve;
      end

      if obj.DEBUGflag
         obj.symFlagLin(obj.nSolve) = symValue;
      else
         obj.symFlagLin(obj.nSolve) = globalsymm;
      end
      if obj.params.iterSinceLastPrecComp == 0
         obj.precCompLin(obj.nSolve) = T_setup;
      else
         obj.precCompLin(obj.nSolve) = 0;
      end
   end

   % Did not converge, if prec not computed for it try again
   if(flag == 1 && obj.params.iterSinceLastPrecComp > 0)
      gresLog().log(3,'Trying to recompute the preconditioner to see if it manages to converge\n');
      obj.params.iterSinceLastPrecComp = 0;
      obj.requestPrecComp = true;
      [x,flag] = obj.SolveLin(A,b,time);
      return;
   end

   % Interesting problem
   if(flag == 1)
      gresLog().log(3,'Iterations since last preconditioner computation %d\n',obj.params.iterSinceLastPrecComp);
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

   gresLog().log(3,'Fallback to matlab due to size or chronos inexistance\n');

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

      diffnorm = norm(A-A','f');
      Anorm = norm(A,'f');
      relNorm = diffnorm/Anorm;

      if relNorm < eps1
         maxval = 0;
         globalsymm = 1;
      else
         maxval = relNorm;
         globalsymm = 0;
      end
      symMat = globalsymm;
      return
   end
   
   % Allocate the stuff
   N = size(A,1);
   symm = ones(sum(1:N),1);
   val = zeros(sum(1:N),1);
   cont = 1;
   
   % Loop over the blocks
   for j = 1:N
      for i = 1:j
         if i == j
            % Diagonal Block
            [symm(cont), val(cont)] = checkSymmetry(A{i,i},eps1);
         elseif ~isempty(A{i,j})
            % Off-Diagonal Block
            diffnorm = norm(A{i,j}-A{j,i}','f');
            Anorm = norm(A{i,j},'f');
            relNorm = diffnorm/Anorm;
            
            
            if relNorm < eps1
                symm(cont) = 1;
                val(cont) = 0;
            else
                symm(cont) = 0;
                val(cont) = relNorm < eps1;
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


function [A] = fixPattern(A)
   N = size(A, 1);
   for j = 1:N
      for i = 1:j
         patt = spones(A{i,j}) - spones(A{j,i}');
         if nnz(patt)
            mask1 = (patt ==  1);  % in A{i,j} but not A{j,i}'
            mask2 = (patt == -1);  % in A{j,i}' but not A{i,j}
            if i ~= j
               A{j,i} = A{j,i} + (A{i,j} .* mask1)' * eps;
               A{i,j} = A{i,j} + (A{j,i} .* mask2')' * eps;

               patt = spones(A{i,j}) - spones(A{j,i}');
               if nnz(patt) ~= 0
                  error('asymmetric pattern found');
               end
            else
               A{i,i} = A{i,i} + (A{i,i} .* mask1)' * eps;
            end
         end
      end
   end
end
