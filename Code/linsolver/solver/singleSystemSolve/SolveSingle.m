function [x, flag, relres, iter, resvec] = SolveSingle(varargin)
% SolveSingle - Overloaded linear solver for single linear systems.
%
% This function provides two overloaded implementations:
%
% -------------------------------------------------------------------------
% SYNTAX 1: Data-Based Implementation (Direct In-Memory Matrices)
% -------------------------------------------------------------------------
%   [x, flag, relres, iter, resvec] = SolveSingle(A, b, TV0Coords, ...
%       ruizFlag, preconditioner, rhsFlag, tol, verb, isCoords, physics, solverType, maxit)
%
%   Inputs:
%     A              - Cell array of matrix blocks (or assembled sparse matrix).
%     b              - Cell array of RHS blocks (or assembled vector).
%     TV0Coords      - Either precomputed test space (e.g. 3N x 6) OR coordinates (N x 3 / N x 2).
%     ruizFlag       - Boolean: Ruiz scaling flag (defaults to true for fixedStress/efemPrec).
%     preconditioner - Preconditioner name ('aAMG', 'aFSAI', 'RACP', 'fixedStress', 'efemPrec', 'none').
%     rhsFlag        - (Optional) 0: input b [default], 1: ones, 2: rand, 3: A*ones.
%     tol            - (Optional) Relative tolerance [default: 1e-6].
%     verb           - (Optional) Verbosity flag [default: true].
%     isCoords       - (Optional) True if TV0Coords contains coordinates [default: auto-detected].
%     physics        - (Optional) Physics string ('displacements', 'pressure', etc.) [default: "displacements"].
%     solverType     - (Optional) 'auto' [default], 'sqmr', 'gmres', or 'direct'.
%     maxit          - (Optional) Maximum iterations [default: 1000].
%
% -------------------------------------------------------------------------
% SYNTAX 2: File-Based Implementation (.mat File Dump)
% -------------------------------------------------------------------------
%   [x, flag, relres, iter, resvec] = SolveSingle(matFile, ...
%       ruizFlag, preconditioner, rhsFlag, tol, verb, isCoords, physics, solverType, maxit)
%
%   Inputs:
%     matFile        - String or char path to .mat file containing A, b, and coordinates.
%     ruizFlag       - (Optional) Overrides ruizFlag from file (if empty, uses file value).
%     preconditioner - (Optional) Overrides preconditioner from file (if empty, uses file value).
%     rhsFlag        - (Optional) 0: input b [default], 1: ones, 2: rand, 3: A*ones.
%     tol            - (Optional) Relative tolerance [default: 1e-6].
%     verb           - (Optional) Verbosity flag [default: true].
%     isCoords       - (Optional) True if coordinates in file [default: true].
%     physics        - (Optional) Physics string [default: "displacements"].
%     solverType     - (Optional) 'auto' [default], 'sqmr', 'gmres', or 'direct'.
%     maxit          - (Optional) Maximum iterations [default: 1000].
%
% Outputs:
%   x      - Solution vector.
%   flag   - Convergence flag (0 = converged, 1 = no convergence).
%   relres - Final relative residual.
%   iter   - Total iteration count.
%   resvec - Residual history vector.
%
% See also SolveSingleFile.

   if nargin < 1
      error('SolveSingle: Not enough input arguments.');
   end

   firstArg = varargin{1};

   % Check if the first argument is a .mat file
   if (ischar(firstArg) || isstring(firstArg)) && ...
      (endsWith(string(firstArg), '.mat', 'IgnoreCase', true) || isfile(firstArg))
      % Dispatch to file-based implementation
      [x, flag, relres, iter, resvec] = solveSingleFromFile(varargin{:});
   else
      % Dispatch to data-based implementation
      [x, flag, relres, iter, resvec] = solveSingleFromData(varargin{:});
   end

end

% =========================================================================
% IMPLEMENTATION 1: File-Based Solver
% =========================================================================
function [x, flag, relres, iter, resvec] = solveSingleFromFile(matFile, ruizFlag, preconditioner, rhsFlag, tol, verb, isCoords, physics, solverType, maxit)

   if nargin < 1 || isempty(matFile)
      error('solveSingleFromFile: MAT file path must be provided.');
   end

   matFile = char(matFile);
   if ~isfile(matFile)
      error('solveSingleFromFile: MAT file "%s" not found.', matFile);
   end

   % Load variables from .mat file
   matData = load(matFile);
   if ~isfield(matData, 'A') || ~isfield(matData, 'b')
      error('solveSingleFromFile: MAT file "%s" must contain variables "A" and "b".', matFile);
   end

   A = matData.A;
   b = matData.b;

   % Extract coordinates / test space
   if isfield(matData, 'coordinates')
      TV0Coords = matData.coordinates;
   elseif isfield(matData, 'TV0Coords')
      TV0Coords = matData.TV0Coords;
   else
      TV0Coords = [];
   end

   % Use metadata from file if optional arguments are omitted
   if (nargin < 2 || isempty(ruizFlag)) && isfield(matData, 'ruizFlag')
      ruizFlag = matData.ruizFlag;
   elseif nargin < 2
      ruizFlag = [];
   end

   if (nargin < 3 || isempty(preconditioner)) && isfield(matData, 'preconditioner') && ~isempty(matData.preconditioner)
      preconditioner = matData.preconditioner;
   elseif nargin < 3
      preconditioner = [];
   end

   if nargin < 4, rhsFlag = 0; end
   if nargin < 5 || isempty(tol), tol = 1e-6; end
   if nargin < 6, verb = true; end
   if nargin < 7 || isempty(isCoords)
      isCoords = ~isempty(TV0Coords);
   end
   if nargin < 8, physics = "displacements"; end
   if nargin < 9, solverType = 'auto'; end
   if nargin < 10, maxit = 1000; end

   if verb
      fprintf('SolveSingle: Loaded system from "%s"\n', matFile);
   end

   % Forward to the data-based solver implementation
   [x, flag, relres, iter, resvec] = solveSingleFromData(A, b, TV0Coords, ...
      ruizFlag, preconditioner, rhsFlag, tol, verb, isCoords, physics, solverType, maxit);

end

% =========================================================================
% IMPLEMENTATION 2: Data-Based Solver
% =========================================================================
function [x, flag, relres, iter, resvec] = solveSingleFromData(A, b, TV0Coords, ruizFlag, preconditioner, rhsFlag, tol, verb, isCoords, physics, solverType, maxit)

   % -------------------------------------------------------------------------
   % 1. Setup paths
   % -------------------------------------------------------------------------
   if exist('gres_root', 'file') == 2
      rootPath = gres_root;
   else
      thisDir = fileparts(mfilename('fullpath'));
      cur = thisDir;
      while ~isempty(cur) && ~isfolder(fullfile(cur, 'ThirdPartyLibs'))
         parent = fileparts(cur);
         if strcmp(parent, cur), break; end
         cur = parent;
      end
      rootPath = cur;
   end
   chronosDir1 = fullfile(rootPath, 'ThirdPartyLibs', 'ChronosLab', 'sources');
   chronosDir2 = fullfile(rootPath, 'ThirdPartyLibs', 'ChronosLab_', 'sources');
   if isfolder(chronosDir1)
      addpath(genpath(chronosDir1));
   elseif isfolder(chronosDir2)
      addpath(genpath(chronosDir2));
   end

   % -------------------------------------------------------------------------
   % 2. Process optional arguments and defaults
   % -------------------------------------------------------------------------
   if nargin < 3
      TV0Coords = [];
   end
   if nargin < 4
      ruizFlag = [];
   end
   if nargin < 5 || isempty(preconditioner)
      preconditioner = 'aAMG';
   end
   if nargin < 6 || isempty(rhsFlag)
      rhsFlag = 0;
   end
   if nargin < 7 || isempty(tol)
      tol = 1e-6;
   end
   if nargin < 8 || isempty(verb)
      verb = true;
   end
   if nargin < 9
      isCoords = [];
   end
   if nargin < 10 || isempty(physics)
      physics = "displacements";
   end
   physics = string(physics);

   if nargin < 11 || isempty(solverType)
      solverType = 'auto';
   end
   if nargin < 12 || isempty(maxit)
      maxit = 1000;
   end

   % Set default ruizFlag if empty
   if isempty(ruizFlag)
      if any(strcmpi(preconditioner, {'fixedstress', 'fixed_stress', 'efemprec', 'efem'}))
         ruizFlag = true;
      else
         ruizFlag = false;
      end
   end

   % -------------------------------------------------------------------------
   % 3. Fix pattern and analyze symmetry
   % -------------------------------------------------------------------------
   [A] = fixPattern(A);
   nsyTol = 100 * eps;
   [globalsymm, maxval, symMat] = checkSymmetry(A, nsyTol);

   if strcmpi(solverType, 'direct') || any(strcmpi(preconditioner, {'direct', 'matlab'}))
      solverType = 'direct';
      ruizFlag = false;
   elseif globalsymm == 0
      if verb
         fprintf('SolveSingle: Matrix is nonsymmetric (max nonsymmetry: %e). Forcing GMRES.\n', maxval);
      end
      solverType = 'gmres';
   elseif strcmpi(solverType, 'auto')
      solverType = 'sqmr';
   else
      solverType = lower(solverType);
   end

   % -------------------------------------------------------------------------
   % 4. Assemble unscaled matrix and construct RHS according to rhsFlag
   % -------------------------------------------------------------------------
   if iscell(A)
      Amat_unscaled = cell2matrix(A);
   else
      Amat_unscaled = A;
   end
   n = size(Amat_unscaled, 1);

   switch rhsFlag
      case 0
         if iscell(b)
            b = -cell2matrix(b);
         end
         if isempty(b)
            error('SolveSingle: Input b is empty with rhsFlag = 0.');
         end
      case 1
         b = ones(n, 1);
      case 2
         b = rand(n, 1);
      case 3
         x_exact = ones(n, 1);
         b = Amat_unscaled * x_exact;
      otherwise
         error('SolveSingle: Unsupported rhsFlag %d. Valid flags are 0, 1, 2, 3.', rhsFlag);
   end
   b_unscaled = b;

   % -------------------------------------------------------------------------
   % 5. Auto-detect isCoords if not explicitly specified
   % -------------------------------------------------------------------------
   if isempty(isCoords)
      if ~isempty(TV0Coords)
         if iscell(TV0Coords)
            isCoords = true;
         elseif isnumeric(TV0Coords)
            % If second dimension <= 3 and number of rows < total system size, it's coordinate data
            if size(TV0Coords, 2) <= 3 && size(TV0Coords, 1) < n
               isCoords = true;
            else
               isCoords = false;
            end
         else
            isCoords = false;
         end
      else
         isCoords = false;
      end
   end

   % -------------------------------------------------------------------------
   % 6. Initialize Ruiz Scaling
   % -------------------------------------------------------------------------
   Ruiz = RuizScaling(verb, 10, 1e-2);
   Ruiz.scalingFlag = ruizFlag;

   % -------------------------------------------------------------------------
   % 7. Build mock problem solver context for preconditioners
   % -------------------------------------------------------------------------
   mockSolver = struct();
   mockSolver.dt = 1;
   mockSolver.simparams.linSolverParams = struct();
   mockSolver.interfaces = {};
   mockSolver.nInterf = 0;

   if isCoords && ~isempty(TV0Coords)
      if iscell(TV0Coords)
         mockSolver.nDom = numel(TV0Coords);
         for d = 1:mockSolver.nDom
            mockSolver.domains(d).grid.coordinates = TV0Coords{d};
         end
      else
         mockSolver.nDom = 1;
         mockSolver.domains(1).grid.coordinates = TV0Coords;
      end
   else
      mockSolver.nDom = 1;
      mockSolver.domains(1).grid.coordinates = zeros(max(1, round(n/3)), 3);
   end

   % -------------------------------------------------------------------------
   % 8. Identify, instantiate, and compute the preconditioner
   % -------------------------------------------------------------------------
   Prec = [];
   if ~strcmpi(solverType, 'direct')
      precKey = lower(char(preconditioner));

      switch precKey
         case {'aamg', 'amg'}
            Prec = aAMG(verb, mockSolver, char(physics));

         case {'afsai', 'fsai'}
            Prec = aFSAI(verb);

         case 'racp'
            Prec = RACP(verb, mockSolver, char(physics));

         case {'fixedstress', 'fixed_stress'}
            Prec = fixedStress(verb, mockSolver);

         case {'efemprec', 'efem'}
            Prec = efemPrec(verb, mockSolver, nsyTol);

         case {'none', 'matlab', 'direct', ''}
            Prec = [];

         otherwise
            error('SolveSingle: Unknown preconditioner "%s".', preconditioner);
      end
   end

   % Compute the preconditioner
   startPrecT = tic;
   if ~isempty(Prec)
      Prec.Ruiz = Ruiz;
      mc = metaclass(Prec);
      m = mc.MethodList(strcmp({mc.MethodList.Name}, 'Compute'));

      if isCoords
         if isa(Prec, 'aAMG')
            if contains(physics, 'pressure', 'IgnoreCase', true) || contains(physics, 'u', 'IgnoreCase', true)
               TV0 = ones(size(A{1,1}, 1), 1);
            else
               TV0 = mk_rbm_3d(TV0Coords);
            end
            A = Prec.Compute(A, symMat, TV0);
         else
            % For RACP, fixedStress, efemPrec: coordinates are in mockSolver
            if ~isempty(m) && ~isempty(m.OutputNames)
               A = Prec.Compute(A, symMat);
            else
               Prec.Compute(A, symMat);
            end
         end
      else
         % Precomputed test space passed in TV0Coords (bypasses test space computation)
         if ~isempty(TV0Coords)
            if ~isempty(m) && ~isempty(m.OutputNames)
               A = Prec.Compute(A, symMat, TV0Coords);
            else
               Prec.Compute(A, symMat, TV0Coords);
            end
         else
            if ~isempty(m) && ~isempty(m.OutputNames)
               A = Prec.Compute(A, symMat);
            else
               Prec.Compute(A, symMat);
            end
         end
      end
   else
      % Preconditioner 'none' / 'matlab'
      if ruizFlag
         A = Ruiz.Compute(A);
      end
   end
   precTime = toc(startPrecT);

   % -------------------------------------------------------------------------
   % 9. Assemble scaled matrix, scale RHS and initial guess
   % -------------------------------------------------------------------------
   if iscell(A)
      Amat = cell2matrix(A);
   else
      Amat = A;
   end

   x0 = zeros(size(b));
   x0 = Ruiz.applyDinv(x0);
   b = Ruiz.applyD(b);

   % -------------------------------------------------------------------------
   % 10. Solve linear system
   % -------------------------------------------------------------------------
   startT = tic;

   if strcmpi(solverType, 'direct') || any(strcmpi(preconditioner, {'matlab', 'direct'}))
      % Direct MATLAB backslash solve
      x = Amat \ b;
      flag = 0;
      relres = norm(Amat*x - b) / norm(b);
      iter = 1;
      resvec = relres;
   else
      % Iterative solve using GMRES or SQMR
      M1 = [];
      M2 = [];
      if ~isempty(Prec)
         M1 = Prec.Apply_L;
         M2 = Prec.Apply_R;
      end

      switch solverType
         case 'gmres'
            restart = min(100, n);
            maxit_outer = ceil(maxit / restart);
            [x, flag, relres, iter1, resvec] = gmres_RIGHT(Amat, b, restart, tol, maxit_outer, M1, M2, x0, verb);
            iter = (iter1(1) - 1) * restart + iter1(2);

         case 'sqmr'
            Afun = @(v) Amat * v;
            [x, flag, relres, iter, resvec] = SQMR(Afun, b, tol, maxit, M1, M2, x0, verb);

         otherwise
            error('SolveSingle: Unknown solverType "%s".', solverType);
      end
   end

   solveTime = toc(startT);

   % Sanity check: imaginary components
   if any(~isreal(x), 'all')
      if verb
         warning('SolveSingle: Complex solution returned. Discarding imaginary part.');
      end
      x = real(x);
   end

   % -------------------------------------------------------------------------
   % 11. De-apply Ruiz scaling and compute residual norms
   % -------------------------------------------------------------------------
   % Compute scaled norms before de-scaling x
   if ruizFlag
      norm_b_scaled = norm(b);
      norm_res_scaled = norm(Amat * x - b);
   else
      norm_b_scaled = [];
      norm_res_scaled = [];
   end

   % De-apply Ruiz scaling to solution
   x = Ruiz.applyD(x);

   % Unscaled norms on original system
   norm_b_unscaled = norm(b_unscaled);
   norm_res_unscaled = norm(Amat_unscaled * x - b_unscaled);

   % Always print the execution recap even if verbosity is false / 0
   printRecap(n, solverType, preconditioner, physics, rhsFlag, ruizFlag, ...
      isCoords, tol, relres, norm_b_unscaled, norm_b_scaled, ...
      norm_res_unscaled, norm_res_scaled, iter, maxit, flag, precTime, solveTime, x);

end

% =========================================================================
% Local Helper Functions
% =========================================================================

function [A] = fixPattern(A)
% Ensure symmetric sparsity pattern across blocks
   if ~iscell(A)
      return;
   end
   N = size(A, 1);
   for j = 1:N
      for i = 1:j
         patt = spones(A{i,j}) - spones(A{j,i}');
         if nnz(patt)
            mask1 = (patt ==  1);
            mask2 = (patt == -1);
            if i ~= j
               A{j,i} = A{j,i} + (A{i,j} .* mask1)' * eps;
               A{i,j} = A{i,j} + (A{j,i} .* mask2')' * eps;
               patt = spones(A{i,j}) - spones(A{j,i}');
               if nnz(patt) ~= 0
                  error('SolveSingle/fixPattern: Asymmetric pattern found.');
               end
            else
               A{i,i} = A{i,i} + (A{i,i} .* mask1)' * eps;
            end
         end
      end
   end
end

function [globalsymm, maxval, symMat] = checkSymmetry(A, eps1)
% Check global and block symmetry
   if ~iscell(A)
      diffnorm = norm(A - A', 'f');
      Anorm = norm(A, 'f');
      if diffnorm == 0 || Anorm == 0
         globalsymm = 1;
         maxval = 0;
      else
         relNorm = diffnorm / Anorm;
         if relNorm < eps1
            maxval = 0;
            globalsymm = 1;
         else
            maxval = relNorm;
            globalsymm = 0;
         end
      end
      symMat = globalsymm;
      return;
   end

   N = size(A, 1);
   symm = ones(sum(1:N), 1);
   val = zeros(sum(1:N), 1);
   cont = 1;

   for j = 1:N
      for i = 1:j
         if i == j
            [symm(cont), val(cont)] = checkSymmetry(A{i,i}, eps1);
         elseif ~isempty(A{i,j})
            diffnorm = norm(A{i,j} - A{j,i}', 'f');
            Anorm = 0.5 * (norm(A{i,j}, 'f') + norm(A{j,i}, 'f'));
            if diffnorm == 0 || Anorm == 0
               symm(cont) = 1;
               val(cont) = 0;
            else
               relNorm = diffnorm / Anorm;
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

   symMat = zeros(N, N);
   symMat(triu(true(N))) = symm;
   symMat = symMat + triu(symMat, 1).';

   globalsymm = min(symm);
   maxval = max(val);
end

function printRecap(n, solverType, preconditioner, physics, rhsFlag, ruizFlag, ...
   isCoords, tolReq, relresSolver, norm_b_unscaled, norm_b_scaled, ...
   norm_res_unscaled, norm_res_scaled, iter, maxit, flag, precTime, solveTime, x)

   rhsDesc = {'Input b (as in SolveLin)', ...
              'Vector of ones: ones(n,1)', ...
              'Random vector: rand(n,1)', ...
              'Manufactured solution: A * ones(n,1)'};
   if rhsFlag >= 0 && rhsFlag <= 3
      rhsStr = rhsDesc{rhsFlag + 1};
   else
      rhsStr = sprintf('Custom (rhsFlag=%d)', rhsFlag);
   end

   if isCoords
      tv0Str = 'Coordinates (RBM constructed / domain coords)';
   else
      tv0Str = 'Precomputed test space (TV0 passed directly)';
   end

   if ruizFlag
      ruizStr = 'Enabled';
   else
      ruizStr = 'Disabled';
   end

   if flag == 0
      statusStr = 'CONVERGED';
   else
      statusStr = 'FAILED TO CONVERGE';
   end

   if strcmpi(solverType, 'direct')
      precStr = 'Bypassed (Direct Solver)';
   elseif isempty(preconditioner)
      precStr = 'None (unpreconditioned)';
   else
      precStr = char(preconditioner);
   end

   fprintf('\n');
   fprintf('======================================================================\n');
   fprintf('                       SOLVESINGLE EXECUTION RECAP                    \n');
   fprintf('======================================================================\n');
   fprintf('  Problem Size (DOFs) : %d\n', n);
   fprintf('  Solver              : %s\n', upper(solverType));
   fprintf('  Preconditioner      : %s\n', precStr);
   fprintf('  Physics Mode        : %s\n', char(physics));
   fprintf('  RHS Type (rhsFlag)  : [%d] %s\n', rhsFlag, rhsStr);
   fprintf('  Ruiz Scaling        : %s\n', ruizStr);
   fprintf('  Test Space / Coords : %s\n', tv0Str);
   fprintf('  Tolerance Requested : %.2e\n', tolReq);
   fprintf('  Solver Rel. Resid.  : %.2e\n', relresSolver);
   if ruizFlag && ~isempty(norm_b_scaled)
      fprintf('  RHS Norm (unscaled) : %.2e\n', norm_b_unscaled);
      fprintf('  RHS Norm (scaled)   : %.2e\n', norm_b_scaled);
      fprintf('  Resid. Norm (unscaled): %.2e (||Ax - b||)\n', norm_res_unscaled);
      fprintf('  Resid. Norm (scaled): %.2e (||A_s x_s - b_s||)\n', norm_res_scaled);
   else
      fprintf('  RHS Norm            : %.2e\n', norm_b_unscaled);
      fprintf('  Residual Norm       : %.2e (||Ax - b||)\n', norm_res_unscaled);
   end
   fprintf('  Iterations Taken    : %d / %d\n', iter, maxit);
   fprintf('  Convergence Status  : %s (flag = %d)\n', statusStr, flag);
   fprintf('  Preconditioner Time : %.4f seconds\n', precTime);
   fprintf('  Linear Solve Time   : %.4f seconds\n', solveTime);
   fprintf('  Total Elapsed Time  : %.4f seconds\n', precTime + solveTime);
   if rhsFlag == 3
      err_exact = norm(x - ones(size(x))) / norm(ones(size(x)));
      fprintf('  Manufactured Error  : ||x - 1|| / ||1|| = %.2e\n', err_exact);
   end
   fprintf('======================================================================\n\n');
end

