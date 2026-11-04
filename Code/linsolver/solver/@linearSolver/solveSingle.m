function [x, stats, flag, relres, iter, resvec] = solveSingle(matFile, varargin)
% solveSingle - Standalone linear solver validation for saved linear systems (.mat).
%
% Syntax:
%   [x, stats] = linearSolver.solveSingle(matFile)
%   [x, stats] = linearSolver.solveSingle(matFile, 'ParameterName', ParameterValue, ...)
%   [x, stats, flag, relres, iter, resvec] = linearSolver.solveSingle(...)
%
% Inputs:
%   matFile - String or char path to .mat file containing linear system data.
%
% Optional Name-Value Parameters (defaults inherit from saved file):
%   'preconditioner' - Preconditioner name ('aAMG', 'fixedStress', 'efemPrec', 'RACP', 'aFSAI', 'none', 'direct')
%   'physics'        - Physics string ('pressure', 'displacements', ["displacements", "pressure"])
%   'ruizFlag'       - Boolean: Ruiz scaling flag
%   'solverType'     - Solver algorithm ('auto', 'sqmr', 'gmres', 'direct')
%   'restart'        - Restart value for GMRES (default: file restart or 100)
%   'tol'            - Relative convergence tolerance (default: file tol or 1e-6)
%   'maxit'          - Maximum iterations (default: file maxit or 1000)
%   'rhsFlag'        - 0: file b [default], 1: ones, 2: rand, 3: A*ones
%   'verb'           - Boolean: verbosity flag [default: true]
%
% Outputs:
%   x      - Solution vector.
%   stats  - Structure containing execution statistics, timings, and residuals.
%   flag   - Convergence flag (0 = converged, 1 = no convergence).
%   relres - Final relative residual.
%   iter   - Iterations count.
%   resvec - Residual history vector.

   if nargin < 1 || isempty(matFile)
      error('linearSolver:solveSingle: matFile must be provided.');
   end

   matFile = char(matFile);
   if ~isfile(matFile)
      error('linearSolver:solveSingle: File "%s" not found.', matFile);
   end

   % -------------------------------------------------------------------------
   % 1. Load variables from .mat file
   % -------------------------------------------------------------------------
   matData = load(matFile);
   if ~isfield(matData, 'A') || ~isfield(matData, 'b')
      error('linearSolver:solveSingle: MAT file "%s" must contain "A" and "b".', matFile);
   end

   A = matData.A;
   b = matData.b;

   if isfield(matData, 'coordinates')
      TV0Coords = matData.coordinates;
   elseif isfield(matData, 'TV0Coords')
      TV0Coords = matData.TV0Coords;
   else
      TV0Coords = [];
   end

   % -------------------------------------------------------------------------
   % 2. Extract simulation defaults from saved file
   % -------------------------------------------------------------------------
   if isfield(matData, 'preconditioner') && ~isempty(matData.preconditioner)
      filePrec = matData.preconditioner;
   else
      if iscell(A) && size(A, 1) > 1
         filePrec = 'fixedStress';
      else
         filePrec = 'aAMG';
      end
   end

   if isfield(matData, 'physics') && ~isempty(matData.physics)
      filePhysics = matData.physics;
   else
      if any(strcmpi(filePrec, {'fixedstress', 'fixed_stress'}))
         filePhysics = ["displacements", "pressure"];
      elseif any(strcmpi(filePrec, {'efemprec', 'efem'}))
         filePhysics = ["displacements", "fractureJump"];
      elseif contains(lower(matFile), {'flow', 'richards', 'pressure', 'flux'})
         filePhysics = "pressure";
      else
         filePhysics = "displacements";
      end
   end

   if isfield(matData, 'ruizFlag') && ~isempty(matData.ruizFlag)
      fileRuiz = matData.ruizFlag;
   else
      if any(strcmpi(filePrec, {'fixedstress', 'fixed_stress', 'efemprec', 'efem'}))
         fileRuiz = true;
      else
         fileRuiz = false;
      end
   end

   if isfield(matData, 'solverType') && ~isempty(matData.solverType)
      fileSolverType = matData.solverType;
   else
      fileSolverType = 'auto';
   end

   if isfield(matData, 'restart') && ~isempty(matData.restart)
      fileRestart = matData.restart;
   else
      fileRestart = 100;
   end

   if isfield(matData, 'tol') && ~isempty(matData.tol)
      fileTol = matData.tol;
   else
      fileTol = 1e-6;
   end

   if isfield(matData, 'maxit') && ~isempty(matData.maxit)
      fileMaxit = matData.maxit;
   else
      fileMaxit = 1000;
   end

   if isfield(matData, 'dt') && ~isempty(matData.dt)
      dt = matData.dt;
   else
      dt = 1.0;
   end

   if isfield(matData, 'R')
      R = matData.R;
   else
      R = [];
   end

   if isfield(matData, 'timeRelaxation') && ~isempty(matData.timeRelaxation)
      timeRelaxation = matData.timeRelaxation;
   else
      timeRelaxation = 0;
   end

   % -------------------------------------------------------------------------
   % 3. Parse input arguments (caller overrides simulation defaults)
   % -------------------------------------------------------------------------
   p = inputParser;
   p.KeepUnmatched = true;
   addParameter(p, 'preconditioner', filePrec);
   addParameter(p, 'physics', filePhysics);
   addParameter(p, 'ruizFlag', fileRuiz);
   addParameter(p, 'solverType', fileSolverType);
   addParameter(p, 'restart', fileRestart);
   addParameter(p, 'tol', fileTol);
   addParameter(p, 'maxit', fileMaxit);
   addParameter(p, 'rhsFlag', 0);
   addParameter(p, 'verb', true);

   parse(p, varargin{:});

   preconditioner = p.Results.preconditioner;
   physics = p.Results.physics;
   ruizFlag = p.Results.ruizFlag;
   solverType = p.Results.solverType;
   restart = p.Results.restart;
   tol = p.Results.tol;
   maxit = p.Results.maxit;
   rhsFlag = p.Results.rhsFlag;
   verb = p.Results.verb;

   % -------------------------------------------------------------------------
   % 4. Setup paths and environment
   % -------------------------------------------------------------------------
   rootPath = '';
   if exist('gres_root', 'file') == 2
      try
         rootPath = gres_root();
      catch
      end
   end
   if isempty(rootPath)
      thisDir = fileparts(mfilename('fullpath'));
      cur = thisDir;
      while ~isempty(cur) && ~isfolder(fullfile(cur, 'ThirdPartyLibs'))
         parent = fileparts(cur);
         if strcmp(parent, cur), break; end
         cur = parent;
      end
      rootPath = cur;
      setappdata(0, 'gres_root', rootPath);
   end
   if isempty(getappdata(0, 'gresLog'))
      setappdata(0, 'gresLog', Logger());
   end
   chronosDir1 = fullfile(rootPath, 'ThirdPartyLibs', 'ChronosLab', 'sources');
   chronosDir2 = fullfile(rootPath, 'ThirdPartyLibs', 'ChronosLab_', 'sources');
   if isfolder(chronosDir1)
      addpath(genpath(chronosDir1));
   elseif isfolder(chronosDir2)
      addpath(genpath(chronosDir2));
   end

   % -------------------------------------------------------------------------
   % 5. Unify matrix representation and analyze symmetry
   % -------------------------------------------------------------------------
   if ~iscell(A)
      A = {A};
   end

   [A] = fixPattern(A);
   nsyTol = 100 * eps;
   [globalsymm, maxval, symMat] = checkSymmetry(A, nsyTol);

   if strcmpi(solverType, 'direct') || any(strcmpi(preconditioner, {'direct', 'matlab'}))
      solverType = 'direct';
      ruizFlag = false;
   elseif globalsymm == 0 && (strcmpi('auto', solverType) || strcmpi('sqmr', solverType))
      if verb
         fprintf('linearSolver:solveSingle: Matrix is nonsymmetric (max nonsymmetry: %e). Forcing GMRES.\n', maxval);
      end
      solverType = 'gmres';
   elseif strcmpi(solverType, 'auto')
      solverType = 'sqmr';
   else
      solverType = lower(solverType);
   end

   % -------------------------------------------------------------------------
   % 6. Construct RHS according to rhsFlag
   % -------------------------------------------------------------------------
   Amat_unscaled = cell2matrix(A);
   n = size(Amat_unscaled, 1);

   switch rhsFlag
      case 0
         if iscell(b)
            b = -cell2matrix(b);
         end
         if isempty(b)
            error('linearSolver:solveSingle: Input b is empty with rhsFlag = 0.');
         end
         b = b / norm(b);
      case 1
         b = ones(n, 1);
      case 2
         b = rand(n, 1);
      case 3
         x_exact = ones(n, 1);
         b = Amat_unscaled * x_exact;
      otherwise
         error('linearSolver:solveSingle: Unsupported rhsFlag %d.', rhsFlag);
   end
   b_unscaled = b;

   % -------------------------------------------------------------------------
   % 7. Initialize Ruiz Scaling
   % -------------------------------------------------------------------------
   Ruiz = RuizScaling(verb, 10, 1e-2);
   Ruiz.scalingFlag = ruizFlag;

   % -------------------------------------------------------------------------
   % 8. Build mock problem solver context for preconditioners
   % -------------------------------------------------------------------------
   mockSolver = struct();
   mockSolver.dt = dt;
   mockSolver.simparams.relTol = tol;
   mockSolver.simparams.linSolverParams = struct();
   mockSolver.interfaces = {};

   if iscell(TV0Coords)
      mockSolver.nDom = numel(TV0Coords);
      for d = 1:mockSolver.nDom
         mockSolver.domains(d).grid.coordinates = TV0Coords{d};
      end
   else
      mockSolver.nDom = 1;
      mockSolver.domains(1).grid.coordinates = TV0Coords;
   end

   if iscell(A) && size(A, 1) > mockSolver.nDom
      mockSolver.nInterf = size(A, 1) - mockSolver.nDom;
   else
      mockSolver.nInterf = 0;
   end

   % Setup dummy physics solver for fixedStress
   dummyBiot = struct();
   dummyBiot.R = R;
   dummyBiot.computeRelaxationMatrix = @() []; % Dummy no-op
   mockSolver.domains(1).getPhysicsSolver = @(name) dummyBiot;

   % -------------------------------------------------------------------------
   % 9. Instantiate and compute Preconditioner
   % -------------------------------------------------------------------------
   physStr = string(physics);
   physname = physStr(1);

   Prec = [];
   if ~strcmpi(solverType, 'direct')
      precKey = lower(char(preconditioner));
      switch precKey
         case {'aamg', 'amg'}
            if any(contains(physics, {'pressure', 'u'}, 'IgnoreCase', true))
               Prec = aAMG(verb, mockSolver, "pressure");
            else
               Prec = aAMG(verb, mockSolver, char(physname));
            end

         case {'afsai', 'fsai'}
            Prec = aFSAI(verb);

         case 'racp'
            Prec = RACP(verb, mockSolver, char(physname));

         case {'fixedstress', 'fixed_stress'}
            Prec = fixedStress(verb, mockSolver);

         case {'efemprec', 'efem'}
            Prec = efemPrec(verb, mockSolver, nsyTol);

         case {'none', 'matlab', 'direct', ''}
            Prec = [];

         otherwise
            error('linearSolver:solveSingle: Unknown preconditioner "%s".', preconditioner);
      end
   end

   % Compute the preconditioner
   startPrecT = tic;
   if ~isempty(Prec)
      Prec.Ruiz = Ruiz;
      mc = metaclass(Prec);
      m = mc.MethodList(strcmp({mc.MethodList.Name}, 'Compute'));

      if isa(Prec, 'aAMG')
         if any(contains(physics, {'pressure', 'u'}, 'IgnoreCase', true))
            TV0 = ones(size(A{1,1}, 1), 1);
         else
            if iscell(TV0Coords)
               TV0 = [];
               for d = 1:numel(TV0Coords)
                  TV0 = [TV0; mk_rbm_3d(TV0Coords{d})];
               end
            else
               TV0 = mk_rbm_3d(TV0Coords);
            end
         end
         A = Prec.Compute(A, symMat, TV0);
      else
         % RACP, fixedStress, efemPrec
         if ~isempty(m) && ~isempty(m.OutputNames)
            A = Prec.Compute(A, symMat);
         else
            Prec.Compute(A, symMat);
         end
      end
   else
      if ruizFlag
         A = Ruiz.Compute(A);
      end
   end
   precTime = toc(startPrecT);

   % Add relaxation matrix computation time if fixedStress was used
   if any(strcmpi(preconditioner, {'fixedstress', 'fixed_stress'})) && timeRelaxation > 0
      precTime = precTime + timeRelaxation;
   end

   % -------------------------------------------------------------------------
   % 10. Assemble scaled matrix, scale RHS and compute scaled tolerance
   % -------------------------------------------------------------------------
   if iscell(A)
      Amat = cell2matrix(A);
   else
      Amat = A;
   end

   x0 = zeros(size(b));
   x0 = Ruiz.applyDinv(x0);
   b = Ruiz.applyD(b);

   % Scale the tolerance as done in SolveLin to guarantee unscaled convergence
   scaledTol = tol;
   if ruizFlag && ~isempty(Ruiz.fullD)
      b_unscaled_norm = norm(b_unscaled);
      b_scaled_norm = norm(b);
      if b_scaled_norm > 0 && b_unscaled_norm > 0
         minD = min(Ruiz.fullD);
         scaleFactor = (minD * b_unscaled_norm) / b_scaled_norm;
         scaledTol = max(tol * min(1, scaleFactor), 100 * eps);
      end
   end

   % -------------------------------------------------------------------------
   % 11. Solve linear system
   % -------------------------------------------------------------------------
   startT = tic;

   if strcmpi(solverType, 'direct') || any(strcmpi(preconditioner, {'matlab', 'direct'}))
      x = Amat \ b;
      flag = 0;
      relres = norm(Amat*x - b) / norm(b);
      iter = 1;
      resvec = relres;
   else
      M1 = [];
      M2 = [];
      if ~isempty(Prec)
         M1 = Prec.Apply_L;
         M2 = Prec.Apply_R;
      end

      switch solverType
         case 'gmres'
            restart_val = min(restart, n);
            maxit_outer = ceil(maxit / restart_val);
            [x, flag, relres, iter1, resvec] = gmres_RIGHT(Amat, b, restart_val, scaledTol, maxit_outer, M1, M2, x0, verb);
            iter = (iter1(1) - 1) * restart_val + iter1(2);

         case 'sqmr'
            Afun = @(v) Amat * v;
            [x, flag, relres, iter, resvec] = SQMR(Afun, b, scaledTol, maxit, M1, M2, x0, verb);

         otherwise
            error('linearSolver:solveSingle: Unknown solverType "%s".', solverType);
      end
   end

   solveTime = toc(startT);

   if any(~isreal(x), 'all')
      if verb
         warning('linearSolver:solveSingle: Complex solution returned. Discarding imaginary part.');
      end
      x = real(x);
   end

   % -------------------------------------------------------------------------
   % 12. De-scale and compute residuals
   % -------------------------------------------------------------------------
   if ruizFlag
      relres_scaled = norm(Amat * x - b) / norm(b);
   else
      relres_scaled = [];
   end

   % De-apply Ruiz scaling to solution
   x = Ruiz.applyD(x);

   % Unscaled relative residual
   relres_unscaled = norm(Amat_unscaled * x - b_unscaled) / norm(b_unscaled);

   % -------------------------------------------------------------------------
   % 13. Populate stats structure
   % -------------------------------------------------------------------------
   stats = struct();
   stats.flag = flag;
   stats.iter = iter;
   stats.maxit = maxit;
   stats.tol = tol;
   stats.scaledTol = scaledTol;
   stats.relres = relres;
   stats.relres_unscaled = relres_unscaled;
   stats.relres_scaled = relres_scaled;
   stats.resvec = resvec;
   stats.solverType = solverType;
   stats.preconditioner = preconditioner;
   stats.physics = physics;
   stats.ruizFlag = ruizFlag;
   stats.precTime = precTime;
   stats.solveTime = solveTime;
   stats.totalTime = precTime + solveTime;

   % Print recap
   printRecap(n, solverType, preconditioner, physics, rhsFlag, ruizFlag, ...
      tol, scaledTol, relres, relres_unscaled, relres_scaled, iter, maxit, flag, ...
      precTime, solveTime, stats.totalTime);

end

% =========================================================================
% Local Helper Functions
% =========================================================================

function [A] = fixPattern(A)
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
                  error('linearSolver:solveSingle/fixPattern: Asymmetric pattern found.');
               end
            else
               A{i,i} = A{i,i} + (A{i,i} .* mask1)' * eps;
            end
         end
      end
   end
end

function [globalsymm, maxval, symMat] = checkSymmetry(A, eps1)
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
   tolReq, scaledTol, relresSolver, relresUnscaled, relresScaled, iter, maxit, flag, ...
   precTime, solveTime, totalTime)

   rhsDesc = {'Input b (from file)', ...
              'Vector of ones: ones(n,1)', ...
              'Random vector: rand(n,1)', ...
              'Manufactured solution: A * ones(n,1)'};
   if rhsFlag >= 0 && rhsFlag <= 3
      rhsStr = rhsDesc{rhsFlag + 1};
   else
      rhsStr = sprintf('Custom (rhsFlag=%d)', rhsFlag);
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
   fprintf('  Physics Mode        : %s\n', strjoin(string(physics), ', '));
   fprintf('  RHS Type (rhsFlag)  : [%d] %s\n', rhsFlag, rhsStr);
   fprintf('  Ruiz Scaling        : %s\n', ruizStr);
   fprintf('  Tolerance Requested : %.2e\n', tolReq);
   if ruizFlag && ~isempty(scaledTol) && scaledTol < tolReq
      fprintf('  Scaled Tolerance    : %.2e (effective solver tol)\n', scaledTol);
   end
   fprintf('  Solver Rel. Resid.  : %.2e\n', relresSolver);
   fprintf('  Unscaled Rel. Resid.: %.2e (||Ax - b|| / ||b||)\n', relresUnscaled);
   if ruizFlag && ~isempty(relresScaled)
      fprintf('  Scaled Rel. Resid.  : %.2e (||A_s x_s - b_s|| / ||b_s||)\n', relresScaled);
   end
   fprintf('  Iterations Taken    : %d / %d\n', iter, maxit);
   fprintf('  Convergence Status  : %s (flag = %d)\n', statusStr, flag);
   fprintf('  Preconditioner Time : %.4f seconds\n', precTime);
   fprintf('  Linear Solve Time   : %.4f seconds\n', solveTime);
   fprintf('  Total Elapsed Time  : %.4f seconds\n', totalTime);
   fprintf('======================================================================\n\n');
end
