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
   
   % Check if the variables have been passed
   if nargin < 5 || isempty(nonlinIter)
      nonlinIter = 1;
   end

   if obj.DEBUGflag
      A
   end

   % Chronos does not exist, continue with matlab default
   if ~obj.ChronosFlag || (getGlobalSize(A) < obj.matlabMaxSize)
      [x,flag] = matlab_solve(obj,A,b,time);
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
   
   % Contact has opened a fracture or something similar so amg does not converge well. 
   % Directly recompute the preconditioner
   if isequal(obj.Prec.phys, 1.1)
      if obj.generalsolver.iterConfig > obj.iterConfigOld && obj.iterConfigOld == 1
         % If the SAM is being used then it is better to just compute it
         % wrt computing the whole preconditioner
         if ~obj.useSAM
            obj.requestPrecComp = true;
         else
            obj.SAM.requestComp = true;
         end
         obj.iterConfigOld = obj.generalsolver.iterConfig;
      elseif obj.generalsolver.iterConfig < obj.iterConfigOld
         obj.iterConfigOld = obj.generalsolver.iterConfig;
      end
   end

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
         T_setup = obj.Prec.getRecomputeSetupTime(T_setup_wall);
      else
         T_setup = T_setup_wall;
      end

      % Update stats
      obj.aTimeComp = obj.aTimeComp + T_setup_wall;
      obj.nComp = obj.nComp + 1;
      obj.whenComputed(length(obj.whenComputed) + 1) = time;
      obj.params.nSolveSinceLastPrecComp = 0;
      gresLog().log(3,'Finished computing the preconditioner\n');

      % Reset the SAM
      obj.SAM.reset();
   else

      % Condense the domains in a 2x2 matrix
      A = obj.Prec.condenseDomains(A);

      % Scale the matrix and the rhs using the precomputed diagonal scaling
      A = obj.Ruiz.scaleMat(A);

      % Update the state of the blocks in certain preconditioners
      obj.Prec.updateStateBlocks(A);

      % Update stats
      obj.params.nSolveSinceLastPrecComp = obj.params.nSolveSinceLastPrecComp + 1;
      T_setup = 0;
   end

   % If the matrix is nonSymmetric then use always GMRES
   if globalsymm == 0
      obj.SolverType = 'gmres';
      gresLog().log(4,'The matrix is nonsymmetric with a maximum nonsymmetry of %e\n',maxval);
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
   obj.x0 = obj.Ruiz.applyDinv(obj.x0);
   b = obj.Ruiz.applyD(b);

   % Compute if necessary the tolerance for the linear solver
   obj.convStrat.computeTol(b, isLinear, nonlinIter, obj.Ruiz, backupb);

   startT = tic;
   obj.params.actualSolver = obj.SolverType;
   switch obj.SolverType
      case 'gmres'

         % Solve the system by GMRES
         [x,flag,obj.params.lastRelres,iter1,obj.params.resvec] = gmres_RIGHT(Amat,b,obj.params.restart,obj.convStrat.Tol,...
                                                            obj.params.maxit/obj.params.restart,...
                                                            obj.SAM.Apply_L,obj.SAM.Apply_R,obj.x0,obj.DEBUGflag);
         obj.params.iter = (iter1(1) - 1) * obj.params.restart + iter1(2);
         
      case 'sqmr'

         % Solve the system by SQMR
         Afun = @(x) Amat*x;
         [x,flag,obj.params.lastRelres,obj.params.iter,obj.params.resvec] = SQMR(Afun,b,obj.convStrat.Tol,obj.params.maxit,...
                                                               obj.SAM.Apply_L,obj.SAM.Apply_R,obj.x0,obj.DEBUGflag);
   end

   % Sanity check for all real solutions
   if any(~isreal(x),'all')
      gresLog().warning(3,'wtf')
      x = real(x);
   end

   % De-apply Ruiz scaling if used
   x = obj.Ruiz.applyD(x);
   Tend = toc(startT);

   % Save statistics for profiling or info in general
   fillStats(obj,Tend,time,symValue,T_setup,obj.SAM.tSetup,obj.convStrat.unscaledTol,obj.convStrat.Tol);

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
      if iscell(backupb)
         backupb_len = length(cell2matrix(backupb));
      else
         backupb_len = length(backupb);
      end
      if backupb_len > 1e5
        error('Chronos could not solve. Fallback to Matlab blocked due to excessive system size')
      end
      [x,~] = matlab_solve(obj,backupA,backupb,time);

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

   % Restore the new starting vector
   obj.x0 = zeros(size(x));
end
















% Function to get and fill the stats
function fillStats(obj,Tend,time,symValue,T_setup,setupTSAM,tolReq,tolScaled)

   obj.aTimeSolve = obj.aTimeSolve + Tend;
   obj.nSolve = obj.nSolve + 1;
   solves = obj.nSolve;
   obj.aIter = obj.aIter + obj.params.iter;
   obj.maxIter = max(obj.maxIter,obj.params.iter);
   obj.precCompLin(solves) = T_setup;
   obj.SAM.CompLin(solves) = setupTSAM;
   obj.iterLin(solves) = obj.params.iter;
   obj.timeLin(solves) = time; 
   obj.solveTLin(solves) = Tend;
   obj.symFlagLin(solves) = symValue;
   obj.tolReqLin(solves) = tolReq;
   obj.scaledTolLin(solves) = tolScaled;
end





% Function to get the global size of the system
function [gSize] = getGlobalSize(A)

   [nBlockRows,nBlockCols] = size(A);
   rowSizes = zeros(nBlockRows,1);
   
   for i = 1:nBlockRows
      if ~isempty(A{i,i})
         rowSizes(i) = size(A{i,i}, 1);
      else
         for j = 1:nBlockCols
            C = A{i,j};
            if ~isempty(C)
               rowSizes(i) = size(C,1);
               break;
            end
         end
      end
   end

   gSize = sum(rowSizes);
end

function [x,flag] = matlab_solve(obj,A,b,time)

   gresLog().log(4,'Fallback to matlab due to size or chronos inexistance\n');

   % Get the full rhs
   if iscell(b)
      b = -cell2matrix(b);
   end

   startT = tic;
   % Solve the system
   A = cell2matrix(A);
   x = A\b;
   Tend = toc(startT);

   % Get the size of the system
   obj.systemSize = max(size(A,1),obj.systemSize);

   symValue = norm(A-A','f')/norm(A,'f');
   obj.params.iter = 0;
   obj.params.actualSolver = 'direct';
   fillStats(obj,Tend,time,symValue,0,0,NaN,NaN);
   obj.Delta_T(obj.nSolve) = 0;
   flag = 0;
end

function [globalsymm,maxval,symMat] = checkSymmetry(A,eps1)
   
   % Base Case: Numeric Matrix
   if ~iscell(A)

      diffnorm = norm(A-A','f');
      Anorm = norm(A,'f');
      
      if diffnorm == 0 || Anorm == 0
         globalsymm = 1;
         maxval = 0;
      else
         relNorm = diffnorm/Anorm;

         if relNorm < eps1
            maxval = 0;
            globalsymm = 1;
         else
            maxval = relNorm;
            globalsymm = 0;
         end
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
         elseif nnz(A{i,j}) > 0 || nnz(A{j,i}) > 0
            % Off-Diagonal Block
            diffnorm = norm(A{i,j}-A{j,i}','f');
            Anorm = 0.5 * (norm(A{i,j},'f') + norm(A{j,i},'f'));
            
            if diffnorm == 0 || Anorm == 0
               symm(cont) = 1;
               val(cont) = 0;
            else
               relNorm = diffnorm/Anorm;

               if relNorm < eps1
                  symm(cont) = 1;
                  val(cont) = 0;
               else
                  symm(cont) = 0;
                  val(cont) = relNorm;
               end
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

function [isLinear] = getIsLinear(generalsolver)
   % Initialize to false to return at the first nonlinear solver in any domain
   isLinear = false;

   % Loop over the different domains
   for i = 1:generalsolver.nDom

      % Get current domain handle
      dom = generalsolver.domains(i);

      % Loop over the various solvers of this domain
      solver = dom.solverNames;
      for j = 1:numel(solver)

         % Call isLinear on the current solver in domain i
         lin = dom.getPhysicsSolver(solver(j)).isLinear();
   
         % Early exit, if one domain has a solver which is nonlinear
         % then all the system is nonlinear
         if lin == false
            return;
         end
      end
   end

   % Loop over the different interfaces
   for i = 1:generalsolver.nInterf

      % Call isLinear on the current interface
      lin = generalsolver.interfaces{i}.isLinear();

      % Early exit, if one interface is nonlinear
      % then all the system is nonlinear
      if lin == false
         return
      end      
   end

   % If reached here all the domains solvers are linear
   isLinear = true;
end

% Ensures block matrix A has non-empty blocks of valid dimensions
% and symmetric nonzero sparsity patterns.
function [A] = fixPattern(A)
   N = size(A, 1);

   % Determine row and column dimensions for each block row/col index.
   % Also flag blocks that are empty [] or 0x0.
   blkRows = zeros(N, 1);
   blkCols = zeros(N, 1);
   emptyMask = false(N, N);

   for j = 1:N
      for i = 1:N
         M = A{i,j};
         if isempty(M)
            emptyMask(i,j) = true;
         else
            [r, c] = size(M);
            if r > 0 && blkRows(i) == 0, blkRows(i) = r; end
            if c > 0 && blkCols(j) == 0, blkCols(j) = c; end
         end
      end
   end

   % For square block partitionings, if an entire row or column had only empty
   % blocks, recover the missing dimension from the symmetric counterpart.
   for i = 1:N
      if blkRows(i) == 0 && blkCols(i) > 0
         blkRows(i) = blkCols(i);
      elseif blkCols(i) == 0 && blkRows(i) > 0
         blkCols(i) = blkRows(i);
      end
   end

   % Replace empty blocks with correctly dimensioned sparse zeros.
   if any(emptyMask, 'all')
      [eRows, eCols] = find(emptyMask);
      for k = 1:numel(eRows)
         i = eRows(k);
         j = eCols(k);
         r = blkRows(i);
         c = blkCols(j);
         if r == 0 || c == 0
            error('fixPattern:cannotDetermineBlockSize', ...
                  'Cannot determine dimensions for empty block (%d,%d).', i, j);
         end
         A{i,j} = sparse(r, c);
      end
   end

   % Symmetrize the sparsity pattern across blocks.
   % If A{i,j} has a nonzero at (r,c) but A{j,i} does not at (c,r), add an
   % eps-weight entry at (c,r) in A{j,i} so the graph structure is symmetric.
   for j = 1:N
      for i = 1:j
         Aij = A{i,j};
         if i == j
            % Diagonal block: check internal symmetry of sparsity pattern
            if nnz(Aij) > 0
               S = spones(Aij);
               St = S';
               if ~isequal(S, St)
                  % Nonzeros in S not present in S'
                  diff = S - St;
                  [r, c] = find(diff == 1);
                  A{i,i} = Aij + sparse(c, r, eps, blkRows(i), blkRows(i));
               end
            end
         else
            % Off-diagonal block pair (i,j) and (j,i)
            Aji = A{j,i};
            nnz_ij = nnz(Aij);
            nnz_ji = nnz(Aji);
            % Fast exit if both blocks are purely zero
            if nnz_ij == 0 && nnz_ji == 0
               continue;
            end

            Sij = spones(Aij);
            Sji_t = spones(Aji)';

            % Fast equality check: compares sparsity structures directly
            if ~isequal(Sij, Sji_t)
               diff = Sij - Sji_t;
               % Entries present in A{i,j} but missing in A{j,i}'
               [r1, c1] = find(diff == 1);
               if ~isempty(r1)
                  A{j,i} = Aji + sparse(c1, r1, eps, blkCols(j), blkRows(i));
               end
               % Entries present in A{j,i}' but missing in A{i,j}
               [r2, c2] = find(diff == -1);
               if ~isempty(r2)
                  A{i,j} = Aij + sparse(r2, c2, eps, blkRows(i), blkCols(j));
               end
            end
         end
      end
   end
end



% Save current linear system matrix, RHS, grid coordinates, and solver metadata to .mat file.
% This enables standalone reproduction and benchmark validation of linear systems via solveSingle.
function saveLinearSystem(obj, A, b, currentSolve, varargin)

   % Determine if matrix export is requested for this solve step
   forcing = ~isempty(varargin) && (isequal(varargin{1}, true) || strcmpi(string(varargin{1}), "true"));
   shouldSave = forcing || strcmpi('all', obj.saveMatType) || ...
                (strcmpi('selected', obj.saveMatType) && ismember(currentSolve, obj.saveMatSolves));
   if ~shouldSave
      return;
   end

   % Extract spatial grid coordinates from domain models
   coordinates = getCoordinates(obj.generalsolver);
   filename = sprintf('linearSystem_solve%d.mat', currentSolve);

   % Extract scaling flag and active preconditioner class
   ruizFlag = ~isempty(obj.Ruiz) && obj.Ruiz.scalingFlag;
   if ~isempty(obj.Prec), preconditioner = class(obj.Prec); else, preconditioner = ''; end

   % For fixed-stress split systems, precompute and record relaxation matrix R.
   % The relaxation matrix (e.g. volumetric fluid coupling in Biot or multiphase poromechanics)
   % is required by fixedStress to decouple mechanics and flow without rebuilding domain matrices.
   timeRelaxation = 0;
   R = [];
   props = struct();
   if isa(obj.Prec, 'fixedStress')
      dom1 = obj.generalsolver.domains(1);
      if ismember("MultiphasePoromechanics", dom1.solverNames)
         biot = dom1.getPhysicsSolver("MultiphasePoromechanics");
      elseif ismember("BiotFullyCoupled", dom1.solverNames)
         biot = dom1.getPhysicsSolver("BiotFullyCoupled");
      else
         error('Configuration not recognized for fixed stress in saveLinearSystem');
      end

      if ~isempty(biot)
         t0_R = tic;
         biot.computeRelaxationMatrix();
         timeRelaxation = toc(t0_R);
         R = biot.R;

         % For multiphase poromechanics, extract fluid saturation and phase densities
         % to enable exact Schur complement flow scaling in standalone solveSingle runs
         if ismember("MultiphasePoromechanics", dom1.solverNames)
            pState = dom1.getState("pressure");
            swState = dom1.getState("saturation");
            cIdx = dom1.dofm.getFieldCells("pressure");
            tpSolver = biot.getFlowSolver();
            sp = pState(cIdx);
            props.sw = swState(cIdx);
            props.rhoW = tpSolver.props.wet.density(sp);
            props.rhoNW = tpSolver.props.nwet.density(sp);
         end
      end
   end

   % Extract physical variable names from domain degree of freedom managers
   if ~isempty(obj.generalsolver.domains)
      physCell = arrayfun(@(x) x.dofm.getVariableNames(), obj.generalsolver.domains, 'UniformOutput', false);
      physics = [physCell{:}];
   elseif ~isempty(obj.Prec) && isa(obj.Prec, 'aAMG')
      if obj.Prec.phys == 0, physics = "pressure"; else, physics = "displacements"; end
   else
      physics = "displacements";
   end

   % Extract solver parameters and time step
   dt = obj.generalsolver.dt;
   solverType = obj.SolverType;
   tol = obj.generalsolver.simparams.relTol;
   maxit = obj.params.maxit;
   if isfield(obj.params, 'restart') && ~isempty(obj.params.restart)
      restart = obj.params.restart;
   else
      restart = 100;
   end

   % Persist system matrix, RHS, and full execution metadata to disk
   save(filename, 'A', 'b', 'coordinates', 'ruizFlag', 'preconditioner', ...
        'physics', 'solverType', 'restart', 'tol', 'maxit', 'dt', 'R', 'timeRelaxation', 'props');
   gresLog().log(3, 'Saved linear system solve %d to %s\n', currentSolve, filename);
end

% Extract spatial grid coordinates from domain objects across single or multiple domains
function coordinates = getCoordinates(generalsolver)
   if numel(generalsolver.domains) == 1
      coordinates = generalsolver.domains(1).grid.coordinates;
   else
      coordinates = arrayfun(@(d) d.grid.coordinates, generalsolver.domains, 'UniformOutput', false);
   end
end
