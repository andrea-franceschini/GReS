function [x,flag] = SolveLin(obj,A,b,time,nonlinIter)
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
% Notes and behavior details:
% - For small problems or when Chronos (external iterative solver) is not
%   available, the method falls back to matlab_solve which uses the direct
%   backslash on the assembled matrix.
% - If the matrix is detected non-symmetric, SolverType is forced to 'gmres'.
% - The preconditioner is computed when requested or when convergence
%   metrics indicate poor performance (high iteration counts or large
%   relative residual). The decision logic is encoded using fields in obj.
% - After a failed solve (flag == 1), the method may attempt to recompute
%   the preconditioner and retry once. If MATLAB direct solve succeeds but
%   Chronos does not, the code falls back onto matlab and continues the simulation.
% - Note that the A passed here might be slightly different than the one
%   passed in the computation of the preconditioner if the two As differ 
%   too much the preconditioner loses effectiveness. Must be recomputed
% - A is passed directly as a cell array, meaning it is already split in 
%   the various blocks (A11,A12,A21,A22 for a single physics single domain
%   with lagrange multipliers)

   if obj.DEBUGflag
      A
   end

   % Chronos does not exist, continue with matlab default
   if ~obj.ChronosFlag || (getGlobalSize(A) < obj.matlabMaxSize)
      [x,flag] = matlab_solve(obj,A,b,time,nonlinIter);
      return
   end

   % Fix the pattern to be symmetric
   [A] = fixPattern(A);
   
   % Check if this step is linear to use maximum resolution needed
   isLinear = getIsLinear(obj.generalsolver);

   % Save a backup for the inputs
   backupA = A;
   backupb = b;

   % Save system to .mat if requested
   saveLinearSystem(obj, backupA, backupb, obj.nSolve + 1);

   % Permute block cell arrays to canonical physics order
   [A, b] = obj.permuteBlocks(A, b);
   
   % Scale the system if required by the preconditioner (e.g. fixedStress rate scaling)
   [A, b] = obj.Prec.scaleSystem(A, b);

   % Check the symmetry of the matrix
   [globalsymm,maxval,symMat] = checkSymmetry(A,obj.nsyTol);

   % Get the full rhs
   if iscell(b)
      b = -cell2matrix(b);
   end

   % Check if the system has changed size and adapt x0 to be of size(b)
   obj.x0 = obj.Prec.checkGrowth(obj,b);
   
   % Contact has had the active set change significantly.
   % Directly recompute the preconditioner
   checkContact(obj);
   
   % Save the solver type
   firstSolver = obj.SolverType;

   % Symmetric matrix and nonsymmetric preconditioner do not work together.
   % There is a risk in finding a spurious nonphysical solution satisfying
   % the residual
   if globalsymm == true && obj.Prec.PrecSym == false
      obj.requestPrecComp = true;
   end
   if globalsymm == true && obj.SAM.PrecSym == false && obj.useSAM
      obj.SAM.requestComp = true;
   end

   % If the matrix is nonSymmetric then use always GMRES
   if globalsymm == 0
      obj.SolverType = 'gmres';
      gresLog().log(4,'The matrix is nonsymmetric with a maximum nonsymmetry of %e\n',maxval);
   end

   % Have the linear solver compute the Preconditioner if necessary
   if(obj.requestPrecComp) 
      gresLog().log(3,'Computing the preconditioner\n');

      time_start = tic;
      % Compute the preconditioner
      A = obj.Prec.Compute(A,symMat);
      T_setup_wall = toc(time_start);

      % If the preconditioner can provide a recomputation-specific setup time
      % (e.g. fixedStress with linear mechanics, where AMGMech is recycled
      % and only the flow part is recomputed), use that for the convergence strategy
      if ismethod(obj.Prec, 'getRecomputeSetupTime')
         T_setup_eff = obj.Prec.getRecomputeSetupTime(T_setup_wall);
      else
         T_setup_eff = T_setup_wall;
      end
      T_setup = T_setup_wall;

      % Update stats
      obj.aTimeComp = obj.aTimeComp + T_setup_wall;
      obj.nComp = obj.nComp + 1;
      obj.whenComputed(length(obj.whenComputed) + 1) = time;
      obj.params.nSolveSinceLastPrecComp = 0;
      gresLog().log(3,'Finished computing the preconditioner\n');

      % Reset the SAM
      obj.SAM.reset();
   else

      % Scale the matrix using the precomputed diagonal scaling
      A = obj.Ruiz.scaleMat(A);

      % Condense the domains and interfaces
      A = obj.Prec.condenseDomains(A);

      % Update the state of the blocks in certain preconditioners
      obj.Prec.updateStateBlocks(A);

      % Update stats
      obj.params.nSolveSinceLastPrecComp = obj.params.nSolveSinceLastPrecComp + 1;
      T_setup = 0;
      T_setup_eff = 0;
   end

   % Convert the matrix to a sparse double if not already like this
   if iscell(A)
      Amat = cell2matrix(A);
   else
      Amat = A;
   end

   % Get the size of the system
   obj.systemSize = size(Amat,1);

   % Get symmetry
   symValue = norm(Amat-Amat','f')/norm(Amat,'f');

   % Store the matrix for the SAM when the preconditioner is being computed
   % anew if SAM is used
   if obj.requestPrecComp
      obj.SAM.getOldMat(Amat,globalsymm);
   end

   % Computes the SAM if deemed necessary and pass the Prec to setup the
   % apply left/right
   obj.SAM.Compute(Amat,globalsymm,obj.Prec);
   
   % Adjust the preconditioner to be of the correct size in the case of
   % growing mesh simulation
   obj.Prec.updateGrowingPrec(Amat);

   % Apply Ruiz scaling if necessary
   b = obj.Ruiz.applyD(b);

   % Compute if necessary the tolerance for the linear solver
   obj.convStrat.computeTol(b, isLinear, nonlinIter, obj.Ruiz, backupb);

   % Solve the linear system with Krylov solvers
   startT = tic;
   obj.params.actualSolver = obj.SolverType;
   [x,flag,obj.params.iter] = krylovSolve(obj,Amat,b);

   % De-apply Ruiz scaling if used
   obj.x0 = x;
   x = obj.Ruiz.applyD(x);

   % Restore solution vector to raw discretizer order
   x = obj.permuteBackBlocks(x);
   Tend = toc(startT);

   % Save statistics for profiling or info in general
   fillStats(obj,Tend,time,symValue,T_setup,T_setup_eff,obj.SAM.tSetup,obj.convStrat.unscaledTol,obj.convStrat.Tol,nonlinIter,obj.params.lastRelres);

   % Did not converge, if the preconditioner was not computed for this matrix
   % recompute it and try again
   if(flag == 1 && obj.params.nSolveSinceLastPrecComp > 0)
      gresLog().log(3,'Trying to recompute the preconditioner to see if it manages to converge\n');
      obj.params.nSolveSinceLastPrecComp = 0;
      obj.requestPrecComp = true;

      [x,flag] = obj.SolveLin(backupA,backupb,time,nonlinIter);
      return;
   end

   % Interesting problem
   if(flag == 1)
      gresLog().log(3,'Number of solves since last preconditioner computation %d\n',obj.params.nSolveSinceLastPrecComp);
      backupb_len = size(cell2matrix(backupb),1);
      if backupb_len > 1e5
        error('Chronos could not solve. Fallback to Matlab blocked due to excessive system size')
      end
      [x,~] = matlab_solve(obj,backupA,backupb,time,nonlinIter);

      % Fall back to matlab to continue without breaking the
      % simulation
      warning('Matlab could solve and chronos did not. Falling back to matlab to continue simulation');
      return;
   end

   % Reset the solver
   obj.SolverType = firstSolver;

   % Check if the preconditioner needs to be recomputed
   obj.convStrat.recomputePrec(obj,Tend);

   % Check if the SAM needs to be (re)computed
   obj.convStrat.recomputeSAM(obj,obj.SAM,Tend);
end
