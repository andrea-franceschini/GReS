classdef hypredrive < precRecycle
   % Adaptive preconditioner reuse strategy translated from
   % hypre-space/hypredrive (PR #173 and PR #177). Decides preconditioner
   % recomputation from rolling solve-time overhead and iteration-ratio
   % stability scores.

   properties (SetAccess = private, GetAccess = public)
      % Match hypredrive/GEOS baseline quirks (bootstrap setup mean, rounded iters)
      hypredriveParity = true

      % Bootstrap and guard parameters
      bootstrapSolves = 3
      minHistoryPoints = 3
      minReuseSolves = 0
      maxReuseSolves = -1
      badDecisionsToRebuild = 2
      maxIterationRatio = 3.0
      maxSolveTimeRatio = -1.0

      % Adaptive scoring parameters
      rebuildThreshold = 1.0
      positiveFloor = 1e-12
      maxObservations = 256

      % Efficiency component (solve_overhead_vs_setup)
      effWeight = 1.0
      effTarget = 1.0
      effScale = 1.0
      amortizationWindow = 10
      effMaxPoints = 5

      % Stability component (iterations RMS ratio to baseline)
      stabWeight = 1.0
      stabTarget = 1.5
      stabScale = 0.5
      stabMaxPoints = 5
   end

   properties (Access = private)
      % Rolling observation buffers for active preconditioner lifetime
      obsIters = []
      obsSolveT = []
      obsSetupT = []

      % Bootstrap baseline state
      baselineValid = false
      baselineIters = 0
      baselineSolveT = 0
      baselineSetupT = 0

      % Consecutive bad-decision hysteresis counter
      badDecisionStreak = 0

      % Per-solve diagnostic history
      stats = struct('time', [], 'iter', [], 'solveT', [], ...
                     'effAgg', [], 'effContrib', [], ...
                     'stabAgg', [], 'stabContrib', [], ...
                     'score', [], 'streak', [], 'status', {{}})
   end

   methods
      function obj = hypredrive(generalsolver)
         % Initialize hypredrive recycling strategy and parse optional hypredriveParity flag
         obj@precRecycle(generalsolver);

         lsp = generalsolver.simparams.linSolverParams;
         if isfield(lsp, 'hypredriveParity')
            obj.hypredriveParity = parseBooleanFlag(lsp.hypredriveParity, 'hypredriveParity');
         end
      end

      function recomputePrec(obj, linsolver, Tend)
         % Evaluate hypredrive adaptive reuse rules and update requestPrecComp
         nsolve = linsolver.nSolve;
         iter = linsolver.params.iter;

         if linsolver.alpha <= 0.
            % Non-positive alpha disables reuse and forces rebuild every solve
            obj.resetState();
            linsolver.Delta_T(nsolve) = 0;
            linsolver.requestPrecComp = true;
            obj.recordStats(linsolver, Tend, NaN, 0, NaN, 0, 0, 'ForceAlpha');
            return
         end

         % Start a new preconditioner lifetime on rebuild; otherwise append reused solve
         if linsolver.requestPrecComp
            linsolver.requestPrecComp = false;
            linsolver.cumTSolveAfterPrec = 0;
            obj.resetState();
            tSetup = linsolver.precCompLin(nsolve);
            obj.recordObservation(iter, Tend, tSetup);
         else
            obj.recordObservation(iter, Tend, 0);
         end

         age = numel(obj.obsIters);
         score = 0;
         effAgg = NaN;
         effContrib = 0;
         stabAgg = NaN;
         stabContrib = 0;

         % Evaluate hard guardrails before adaptive scoring
         if obj.maxReuseSolves >= 0 && age >= obj.maxReuseSolves
            linsolver.requestPrecComp = true;
            status = 'GuardMaxAge';
         elseif obj.baselineValid && obj.maxIterationRatio > 0 && ...
                obj.baselineIters > 0 && (iter / obj.baselineIters) >= obj.maxIterationRatio
            linsolver.requestPrecComp = true;
            status = 'GuardIterRatio';
         elseif obj.baselineValid && obj.maxSolveTimeRatio > 0 && ...
                obj.baselineSolveT > 0 && (Tend / obj.baselineSolveT) >= obj.maxSolveTimeRatio
            linsolver.requestPrecComp = true;
            status = 'GuardTimeRatio';
         % Defer scoring during bootstrap, minimum reuse age, or insufficient history
         elseif ~obj.baselineValid
            obj.badDecisionStreak = 0;
            status = 'Bootstrap';
         elseif obj.minReuseSolves > 0 && age < obj.minReuseSolves
            obj.badDecisionStreak = 0;
            status = 'GuardMinAge';
         elseif age < obj.minHistoryPoints
            obj.badDecisionStreak = 0;
            status = 'InsufficientHist';
         else
            % Mean clamped solve-time overhead relative to amortized setup budget
            effIdx = max(1, age - obj.effMaxPoints + 1):age;
            budget = max(obj.baselineSetupT / obj.amortizationWindow, obj.positiveFloor);
            effSamples = max(obj.obsSolveT(effIdx) - obj.baselineSolveT, 0) / budget;
            effAgg = mean(effSamples);
            effContrib = obj.effWeight * max(0, (effAgg - obj.effTarget) / obj.effScale);

            % RMS of iteration ratio relative to baseline
            stabIdx = max(1, age - obj.stabMaxPoints + 1):age;
            % Parity mode rounds baselineIters to integer to match hypredrive C code
            if obj.hypredriveParity
               baseIter = max(round(obj.baselineIters), obj.positiveFloor);
            else
               baseIter = max(obj.baselineIters, obj.positiveFloor);
            end
            stabSamples = obj.obsIters(stabIdx) / baseIter;
            stabAgg = sqrt(mean(stabSamples .^ 2));
            stabContrib = obj.stabWeight * max(0, (stabAgg - obj.stabTarget) / obj.stabScale);

            % Require badDecisionsToRebuild consecutive bad scores before requesting rebuild
            score = effContrib + stabContrib;
            isBad = (score >= obj.rebuildThreshold);
            if isBad
               obj.badDecisionStreak = obj.badDecisionStreak + 1;
            else
               obj.badDecisionStreak = 0;
            end

            if isBad && obj.badDecisionStreak >= obj.badDecisionsToRebuild
               linsolver.requestPrecComp = true;
               status = 'Rebuild';
            elseif isBad
               status = 'Hold';
            else
               status = 'Good';
            end
         end

         % Store composite score in Delta_T for unified linearSolver diagnostics
         linsolver.Delta_T(nsolve) = score;
         obj.recordStats(linsolver, Tend, effAgg, effContrib, stabAgg, stabContrib, score, status);
      end

      function recomputeSAM(~, linsolver, ~, ~)
         % Reject SAM explicitly since hypredrive does not model SAM recycling
         if linsolver.useSAM
            error('hypredrive:SAMNotSupported', ...
               'SAM is not supported with the hypredrive recycle strategy.');
         end
      end

      function printStats(obj)
         % Print per-solve hypredrive scoring table and status distribution
         if isempty(obj.stats.time)
            return
         end

         lineLen = 109;
         sepLine = [repmat('-', 1, lineLen), '\n'];
         fprintf(['\n', sepLine]);
         fprintf('Hypredrive Adaptive Reuse Statistics (hypredriveParity=%d)\n', obj.hypredriveParity);
         fprintf(sepLine);
         fprintf('| %8s | %6s | %5s | %8s | %8s | %8s | %8s | %6s | %-15s |\n', ...
            'PhysTime', 'Sol N.', 'Iter', 'SolveT', 'Eff', 'Stab', 'Score', 'Streak', 'Status');
         fprintf(sepLine);

         n = length(obj.stats.time);
         for i = 1:n
            effStr = formatOptionalFloat(obj.stats.effAgg(i));
            stabStr = formatOptionalFloat(obj.stats.stabAgg(i));
            fprintf('| %.2e | %6d | %5d | %.2e | %8s | %8s | %.2e | %6d | %-15s |\n', ...
               obj.stats.time(i), i, obj.stats.iter(i), obj.stats.solveT(i), ...
               effStr, stabStr, obj.stats.score(i), obj.stats.streak(i), obj.stats.status{i});
         end
         fprintf(sepLine);

         statusList = obj.stats.status;
         uniqueStatuses = unique(statusList);
         parts = cell(1, numel(uniqueStatuses));
         for s = 1:numel(uniqueStatuses)
            st = uniqueStatuses{s};
            parts{s} = sprintf('%s=%d', st, sum(strcmp(statusList, st)));
         end
         fprintf('  Status counts: %s\n', strjoin(parts, ', '));
         fprintf(sepLine);
         fprintf('\n');
      end
   end

   methods (Access = private)
      function resetState(obj)
         % Clear rolling observations, bootstrap baseline, and streak for a new prec
         obj.obsIters = [];
         obj.obsSolveT = [];
         obj.obsSetupT = [];
         obj.baselineValid = false;
         obj.baselineIters = 0;
         obj.baselineSolveT = 0;
         obj.baselineSetupT = 0;
         obj.badDecisionStreak = 0;
      end

      function recordObservation(obj, iter, solveT, setupT)
         % Append solve metrics to rolling window and lock baseline after bootstrapSolves
         obj.obsIters(end + 1) = iter;
         obj.obsSolveT(end + 1) = solveT;
         obj.obsSetupT(end + 1) = setupT;

         if numel(obj.obsIters) > obj.maxObservations
            obj.obsIters(1) = [];
            obj.obsSolveT(1) = [];
            obj.obsSetupT(1) = [];
         end

         if ~obj.baselineValid && numel(obj.obsIters) >= obj.bootstrapSolves
            obj.baselineIters = mean(obj.obsIters(1:obj.bootstrapSolves));
            obj.baselineSolveT = mean(obj.obsSolveT(1:obj.bootstrapSolves));
            % Parity mode averages setupT over bootstrap window (including zeros); non-parity uses initial setupT
            if obj.hypredriveParity
               obj.baselineSetupT = mean(obj.obsSetupT(1:obj.bootstrapSolves));
            else
               obj.baselineSetupT = obj.obsSetupT(1);
            end
            obj.baselineValid = true;
         end
      end

      function recordStats(obj, linsolver, solveT, effAgg, effContrib, stabAgg, stabContrib, score, status)
         % Record per-solve scores, hysteresis streak, and decision status for printStats
         idx = linsolver.nSolve;
         obj.stats.time(idx)        = linsolver.timeLin(idx);
         obj.stats.iter(idx)        = linsolver.params.iter;
         obj.stats.solveT(idx)      = solveT;
         obj.stats.effAgg(idx)      = effAgg;
         obj.stats.effContrib(idx)  = effContrib;
         obj.stats.stabAgg(idx)     = stabAgg;
         obj.stats.stabContrib(idx) = stabContrib;
         obj.stats.score(idx)       = score;
         obj.stats.streak(idx)      = obj.badDecisionStreak;
         obj.stats.status{idx}      = status;
      end
   end
end

function flag = parseBooleanFlag(val, name)
   % Parse logical, 0/1 numeric, or "true"/"false"/"1"/"0" string flag strictly
   if islogical(val) && isscalar(val)
      flag = val;
   elseif isnumeric(val) && isscalar(val) && (val == 0 || val == 1)
      flag = logical(val);
   elseif (ischar(val) || isstring(val)) && isscalar(string(val))
      cleanVal = strtrim(string(val));
      if strcmpi(cleanVal, "true") || cleanVal == "1"
         flag = true;
      elseif strcmpi(cleanVal, "false") || cleanVal == "0"
         flag = false;
      else
         error('hypredrive:InvalidInput', ...
            'Expected "%s" to be true/false or 1/0, but received: %s', name, cleanVal);
      end
   else
      error('hypredrive:InvalidInputType', ...
         'Invalid type or dimension for "%s". Expected scalar boolean, 0/1 numeric, or valid string/char.', name);
   end
end

function str = formatOptionalFloat(val)
   % Format scalar float or placeholder dashes when NaN
   if isnan(val)
      str = '     ---';
   else
      str = sprintf('%8.2e', val);
   end
end
