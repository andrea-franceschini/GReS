function [x, stats, flag, relres, iter, resvec] = solveSingle(matFile, varargin)
% solveSingle - Standalone linear solver execution for saved linear systems (.mat).
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
%   'preconditioner' - Preconditioner name ('aAMG', 'fixedStress', 'CPR', 'efemPrec', 'RACP', 'aFSAI', 'none', 'direct')
%   'physics'        - Physics string ('pressure', 'displacements', ["displacements", "pressure"], etc.)
%   'ruizFlag'       - Boolean: Ruiz scaling flag
%   'solverType'     - Solver algorithm ('auto', 'sqmr', 'gmres', 'direct', 'matlab')
%   'useMatlab'      - Boolean: force MATLAB direct backslash solve [default: false]
%   'useIterative'   - Boolean: force iterative solver even for small matrices [default: true]
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
%   relres - Final relative residual reported by solver.
%   iter   - Iteration count.
%   resvec - Residual history vector.

   if nargin < 1 || isempty(matFile)
      error('linearSolver:solveSingle: matFile must be provided.');
   end

   matFile = char(matFile);
   if ~isfile(matFile)
      error('linearSolver:solveSingle: File "%s" not found.', matFile);
   end

   % Load linear system data from mat file
   matData = load(matFile);
   if ~isfield(matData, 'A') || ~isfield(matData, 'b')
      error('linearSolver:solveSingle: MAT file "%s" must contain "A" and "b".', matFile);
   end

   A = matData.A;
   b = matData.b;

   % Parse inputs with simulation defaults inheriting directly from file metadata
   p = inputParser;
   p.KeepUnmatched = true;
   addParameter(p, 'preconditioner', getFieldOr(matData, 'preconditioner', 'aAMG'));
   addParameter(p, 'physics',        getFieldOr(matData, 'physics', "displacements"));
   addParameter(p, 'ruizFlag',       getFieldOr(matData, 'ruizFlag', false));
   addParameter(p, 'solverType',     getFieldOr(matData, 'solverType', 'auto'));
   addParameter(p, 'restart',        getFieldOr(matData, 'restart', 100));
   addParameter(p, 'tol',            getFieldOr(matData, 'tol', 1e-6));
   addParameter(p, 'maxit',          getFieldOr(matData, 'maxit', 1000));
   addParameter(p, 'useMatlab',      false);
   addParameter(p, 'useIterative',   true);
   addParameter(p, 'rhsFlag',        0);
   addParameter(p, 'verb',           true);

   parse(p, varargin{:});
   opt = p.Results;

   coordinates = getFieldOr(matData, 'coordinates', []);
   dt          = getFieldOr(matData, 'dt', 1.0);
   R           = getFieldOr(matData, 'R', []);

   % Ensure block matrix representation and construct RHS vector
   if ~iscell(A), A = {A}; end
   Amat_unscaled = cell2matrix(A);
   n = size(Amat_unscaled, 1);
   [b_in, b_unscaled] = constructRhs(b, Amat_unscaled, A, opt.rhsFlag);

   % Extract fluid properties for multiphase fixed-stress systems.
   % For legacy matrices lacking saved fluid properties, fall back to default
   % wetting saturation sw = 0.5 and unit phase densities.
   ncFlow = size(A{min(2, size(A, 1)), min(2, size(A, 2))}, 1);
   props = struct('sw', 0.5 * ones(ncFlow, 1), 'rhoW', ones(ncFlow, 1), 'rhoNW', ones(ncFlow, 1));
   if isfield(matData, 'props')
      props = matData.props;
   end

   % Determine whether to force MATLAB direct backslash solve.
   % For direct solves, Ruiz diagonal equilibration and iterative preconditioners
   % are bypassed because backslash handles factorization directly.
   isDirect = opt.useMatlab || any(strcmpi(opt.solverType, {'direct', 'matlab'})) || ...
              any(strcmpi(opt.preconditioner, {'direct', 'matlab'}));

   linParams = struct();
   if isDirect
      linParams.useMatlab = 1;
      opt.ruizFlag = false;
   else
      linParams.useIterative = 1;
      linParams.scalingFlag = opt.ruizFlag;
      linParams.tol = opt.tol;
      linParams.maxit = opt.maxit;
      linParams.restart = opt.restart;
   end

   % Construct mock simulation context providing domain DOF managers and physics solvers.
   % This allows linearSolver to automatically choose physics-aware preconditioners (CPR,
   % fixedStress, RACP) without duplicating assembly or simulation harness code.
   mockSolver = buildMockSolver(coordinates, opt.physics, A, dt, opt.tol, linParams, R, props);
   lin = linearSolver(mockSolver, opt.physics);

   % Apply explicit user solver and preconditioner overrides
   if ~isDirect
      if ~strcmpi(opt.solverType, 'auto')
         lin.SolverType = lower(opt.solverType);
      end
      if any(strcmpi(opt.preconditioner, {'afsai', 'fsai'}))
         lin.Prec = aFSAI(opt.verb);
         lin.ChronosFlag = true;
      elseif any(strcmpi(opt.preconditioner, {'none', 'unpreconditioned'}))
         lin.Prec = [];
      end
   end

   % Execute linear solver through SolveLin
   [x, flag] = lin.SolveLin(A, b_in, 0, 1);

   % Package statistics and calculate unscaled and true residual norms
   [stats, relres, iter, resvec] = packageStats(lin, x, flag, Amat_unscaled, b_unscaled, ...
      opt.tol, opt.maxit, opt.physics, opt.ruizFlag, isDirect);

   % Print summary table when requested
   if opt.verb
      printRecap(n, stats.solverType, stats.preconditioner, opt.physics, opt.rhsFlag, ...
         opt.ruizFlag, opt.tol, stats.scaledTol, relres, stats.trueResNorm, ...
         stats.relres_unscaled, iter, opt.maxit, flag, stats.precTime, stats.solveTime, stats.totalTime);
   end

end

% =========================================================================
% Local Helper Functions
% =========================================================================

% Helper function to extract a struct field or return a default value
function val = getFieldOr(s, fieldName, defaultVal)
   if isfield(s, fieldName) && ~isempty(s.(fieldName))
      val = s.(fieldName);
   else
      val = defaultVal;
   end
end

% Construct right-hand side vector based on the requested rhsFlag
function [b_in, b_unscaled] = constructRhs(b_file, Amat, A, rhsFlag)
   n = size(Amat, 1);
   switch rhsFlag
      case 0
         % File RHS: if stored as cell array, convert to flat unscaled vector
         if iscell(b_file)
            b_unscaled = -cell2matrix(b_file);
            b_in = b_file;
         else
            b_unscaled = b_file;
            b_in = b_file;
         end
         if isempty(b_unscaled)
            error('linearSolver:solveSingle: Input b is empty with rhsFlag = 0.');
         end
      case 1
         % Uniform unit vector
         b_unscaled = ones(n, 1);
         b_in = partitionRhs(b_unscaled, A);
      case 2
         % Random vector
         b_unscaled = rand(n, 1);
         b_in = partitionRhs(b_unscaled, A);
      case 3
         % Manufactured solution with exact x = ones(n, 1)
         b_unscaled = Amat * ones(n, 1);
         b_in = partitionRhs(b_unscaled, A);
      otherwise
         error('linearSolver:solveSingle: Unsupported rhsFlag %d.', rhsFlag);
   end
end

% Partition flat RHS into a cell array matching block rows of A when multi-physics
function b_in = partitionRhs(b_unscaled, A)
   if iscell(A) && size(A, 1) > 1
      rowSizes = getBlockRowSizes(A);
      % SolveLin internally negates cell RHS arrays (assuming -residual convention),
      % so we negate here to preserve positive RHS for custom manufactured vectors
      b_in = mat2cell(-b_unscaled, rowSizes, 1);
   else
      b_in = b_unscaled;
   end
end

% Compute number of rows per block row, safely ignoring empty off-diagonal blocks
function rowSizes = getBlockRowSizes(A)
   nRows = size(A, 1);
   rowSizes = zeros(nRows, 1);
   for i = 1:nRows
      for j = 1:size(A, 2)
         if ~isempty(A{i, j})
            rowSizes(i) = size(A{i, j}, 1);
            break;
         end
      end
   end
end

% Build a lightweight simulation harness mimicking generalsolver.
% This provides degrees of freedom, coordinate mappings, and physics solver
% handles needed by coupled block preconditioners (fixedStress, CPR, RACP)
% without requiring the full simulation problem or mesh classes.
function mockSolver = buildMockSolver(coordinates, physics, A, dt, tol, linParams, R, props)

   % Normalize coordinates into cell array to treat single-domain and
   % multi-domain problems uniformly without nested cell checks
   if ~iscell(coordinates), coordinates = {coordinates}; end

   % Derive domain and interface counts from block matrix and coordinates.
   % Any block rows beyond the number of spatial domains correspond to interface
   % or Lagrange multiplier constraints (e.g., mortar contact interfaces).
   nDom = numel(coordinates);
   nInt = max(0, size(A, 1) - nDom);

   % Determine variable assignment per domain.
   % In multi-domain single physics (e.g. ["displacements", "displacements"]),
   % each domain holds 1 variable. In single-domain coupled physics
   % (e.g. ["displacements", "pressure"]), domain 1 holds all coupled variables.
   physStr = string(physics);
   if nDom > 1 && isscalar(unique(physStr))
      nVarsPerDom = 1;
      domPhysics = unique(physStr);
   else
      nVarsPerDom = numel(physStr);
      domPhysics = physStr;
   end

   % Construct mock domains with DOF manager interfaces and dummy solvers
   ncFlow = size(A{min(2, size(A, 1)), min(2, size(A, 2))}, 1);
   dom = struct();
   for d = 1:nDom
      coords_d = coordinates{d};
      nc = max(size(coords_d, 1), size(A{d, d}, 1));

      dom(d).grid.coordinates = coords_d;
      dom(d).solverNames = "GenericSolver";
      dom(d).dofm.getNumberOfVariables = @() nVarsPerDom;
      dom(d).dofm.getVariableNames = @() domPhysics;
      dom(d).dofm.getFieldCells = @(name) 1:ncFlow;
      dom(d).dofm.getNumbDoF = @(name) nc;
      dom(d).getState = @(name) getMockState(name, props.sw, ncFlow, nc);

      % Provide physics solver handles for coupled mechanics and flow relaxation
      dummyMech = struct('isLinear', @() true);
      dummyFlow = struct('props', struct('wet', struct('density', @(p) props.rhoW), ...
                                         'nwet', struct('density', @(p) props.rhoNW)));
      dummySolver = struct('R', R, 'computeRelaxationMatrix', @() [], 'isLinear', @() true, ...
                           'getMechSolver', @() dummyMech, 'getFlowSolver', @() dummyFlow);
      dom(d).getPhysicsSolver = @(name) dummySolver;
   end

   % Initialize linear interface placeholders for mortar contact blocks
   interfaces = repmat({struct('isLinear', @() true)}, nInt, 1);

   % Assemble mock solver root structure
   mockSolver = struct();
   mockSolver.dt = dt;
   mockSolver.simparams.relTol = tol;
   mockSolver.simparams.absTol = 1e-12;
   mockSolver.simparams.linSolverParams = linParams;
   mockSolver.nDom = nDom;
   mockSolver.nInterf = nInt;
   mockSolver.interfaces = interfaces;
   mockSolver.domains = dom;
   mockSolver.iterConfig = 1;
   mockSolver.iterConfigOld = 1;
end

% Collect solve statistics, timings, and calculate true and relative residuals
function [stats, relres, iter, resvec] = packageStats(lin, x, flag, Amat_unscaled, b_unscaled, ...
   tol, maxit, physics, ruizFlag, isDirect)

   % Compute unscaled true algebraic residual norm: ||A*x - b||
   norm_b = norm(b_unscaled);
   trueResNorm = norm(Amat_unscaled * x - b_unscaled);
   relres_unscaled = trueResNorm / max(norm_b, eps);

   % Retrieve solver-reported relative residual, iteration count, and residual history
   relres = getFieldOr(lin.params, 'lastRelres', relres_unscaled);
   iter   = getFieldOr(lin.params, 'iter', 0);
   resvec = getFieldOr(lin.params, 'resvec', []);

   % Extract setup and solve timings from solver history records
   precTime  = getLastOr(lin.precCompLin, 0);
   solveTime = getLastOr(lin.solveTLin, 0);
   scaledTol = getLastOr(lin.scaledTolLin, []);

   % Identify active solver and preconditioner labels for reporting.
   % For direct solves, backslash handles factorization directly without preconditioner.
   % For iterative solves, actualSolver reflects the runtime algorithm (SQMR vs GMRES).
   if isDirect
      actualSolver = 'DIRECT (Matlab)';
      actualPrec   = 'Bypassed (Direct Solver)';
   else
      actualSolver = upper(getFieldOr(lin.params, 'actualSolver', lin.SolverType));
      if isempty(lin.Prec)
         actualPrec = 'None (unpreconditioned)';
      else
         actualPrec = class(lin.Prec);
      end
   end

   % Package comprehensive solve statistics struct
   stats = struct();
   stats.flag = flag;
   stats.iter = iter;
   stats.maxit = maxit;
   stats.tol = tol;
   stats.scaledTol = scaledTol;
   stats.relres = relres;
   stats.trueResNorm = trueResNorm;
   stats.relres_unscaled = relres_unscaled;
   stats.resvec = resvec;
   stats.solverType = actualSolver;
   stats.preconditioner = actualPrec;
   stats.physics = physics;
   stats.ruizFlag = ruizFlag;
   stats.precTime = precTime;
   stats.solveTime = solveTime;
   stats.totalTime = precTime + solveTime;
end

% Return last element of an array if non-empty, otherwise return default value
function val = getLastOr(arr, defaultVal)
   if ~isempty(arr)
      val = arr(end);
   else
      val = defaultVal;
   end
end

% Print formatted summary of the solve execution to console
function printRecap(n, solverType, precStr, physics, rhsFlag, ruizFlag, ...
   tolReq, scaledTol, relresSolver, trueResNorm, relresUnscaled, iter, maxit, flag, ...
   precTime, solveTime, totalTime)

   % Map RHS numeric flag to descriptive text
   rhsDesc = {'Input b (from file)', ...
              'Vector of ones: ones(n,1)', ...
              'Random vector: rand(n,1)', ...
              'Manufactured solution: A * ones(n,1)'};
   if rhsFlag >= 0 && rhsFlag <= 3
      rhsStr = rhsDesc{rhsFlag + 1};
   else
      rhsStr = sprintf('Custom (rhsFlag=%d)', rhsFlag);
   end

   % Format Ruiz status and convergence status strings
   if ruizFlag, ruizStr = 'Enabled'; else, ruizStr = 'Disabled'; end
   if flag == 0, statusStr = 'CONVERGED'; else, statusStr = 'FAILED TO CONVERGE'; end

   % Print minimal, clean summary lines without banners
   fprintf('\nsolveSingle summary:\n');
   fprintf('  Problem Size        : %d DOFs\n', n);
   fprintf('  Solver              : %s\n', char(solverType));
   fprintf('  Preconditioner      : %s\n', char(precStr));
   fprintf('  Physics Mode        : %s\n', strjoin(string(physics), ', '));
   fprintf('  RHS Type            : [%d] %s\n', rhsFlag, rhsStr);
   fprintf('  Ruiz Scaling        : %s\n', ruizStr);
   fprintf('  Tolerance Requested : %.2e\n', tolReq);
   if ruizFlag && ~isempty(scaledTol) && scaledTol < tolReq
      fprintf('  Scaled Tolerance    : %.2e\n', scaledTol);
   end
   fprintf('  Solver Rel. Resid.  : %.2e\n', relresSolver);
   fprintf('  True Residual Norm  : %.2e (||Ax - b||)\n', trueResNorm);
   fprintf('  Unscaled Rel. Resid.: %.2e (||Ax - b|| / ||b||)\n', relresUnscaled);
   fprintf('  Iterations Taken    : %d / %d\n', iter, maxit);
   fprintf('  Convergence Status  : %s (flag = %d)\n', statusStr, flag);
   fprintf('  Preconditioner Time : %.4f s\n', precTime);
   fprintf('  Linear Solve Time   : %.4f s\n', solveTime);
   fprintf('  Total Elapsed Time  : %.4f s\n\n', totalTime);
end

% Return mock state vector by physical variable name
function state = getMockState(name, sw, ncFlow, nc)
   if strcmpi(name, "saturation")
      state = sw;
   elseif strcmpi(name, "pressure")
      state = zeros(ncFlow, 1);
   else
      state = zeros(nc, 1);
   end
end
