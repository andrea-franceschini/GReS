classdef EisenstatWalker < tolStrategy
   % Eisenstat-Walker Choice 2 forcing term strategy (Eisenstat & Walker, 1996).

   properties (SetAccess = private, GetAccess = public)
      % Eisenstat-Walker parameters (Choice 2)
      gamma = 0.9
      alpha = 2.0
      maxEtak = 0.01
      minEtak
      oldEtak
      eta0 = 0.005
      normFxold
   end

   properties (Access = private)
      % Retry detection
      lastNonlinIter = -1
      lastPhysNorm = -1
      lastRawEtak = NaN
      lastRatio = NaN

      % Step info for current solve
      curStep = struct('ratio', NaN, 'rawEta', NaN, 'eta', NaN, ...
                       'tolReq', NaN, 'nonlinIter', 0, 'status', '')

      % Historical statistics table
      stats = struct('time', [], 'nonlinIter', [], 'ratio', [], ...
                     'rawEta', [], 'eta', [], 'tolReq', [], ...
                     'tolReached', [], 'iter', [], 'itPerDec', [], ...
                     'deltaT', [], 'status', {{}})
   end

   methods
      function obj = EisenstatWalker(generalsolver)
         % Set floor minEtak = linearTol and read optional initial forcing term eta0
         obj@tolStrategy(generalsolver);
         obj.minEtak = obj.linearTol;

         lsp = generalsolver.simparams.linSolverParams;
         if isfield(lsp, 'eta0')
            obj.eta0 = str2double(string(lsp.eta0));
         end
      end

      function computeTol(obj, b, islinear, nonlinIter, Ruiz, backupb)
         % Compute Eisenstat-Walker Choice 2 forcing term eta_k for the current Newton step
         physNorm = obj.getPhysNorm(backupb);
         % Detect Krylov retry on the same Newton residual (same iter and RHS norm)
         isRetry = (nonlinIter == obj.lastNonlinIter && abs(physNorm - obj.lastPhysNorm) <= 1e-12 * max(1, physNorm));

         if islinear
            % Solve to floor tolerance with Ruiz adjustment
            targetTol = obj.minEtak;
            rawEtak = obj.minEtak;
            ratio = NaN;
            status = 'Linear';

            obj.applyLinearRuizTol(targetTol, b, Ruiz, backupb, physNorm);
         elseif isRetry
            % Preserve tolerance on retry without altering history
            targetTol = obj.unscaledTol;
            rawEtak = obj.lastRawEtak;
            ratio = obj.lastRatio;
            status = 'Retry';

            obj.unscaledTol = targetTol;
            obj.Tol = targetTol;
         else
            if nonlinIter ~= 1
               % Compute Choice 2 forcing term eta_k = gamma * (||F_k|| / ||F_{k-1}||)^alpha
               ratio = physNorm / obj.normFxold;
               rawEtak = obj.gamma * (ratio)^obj.alpha;
               etak = rawEtak;
               status = 'Choice2';

               % Prevent etak from dropping too rapidly (Choice 2 safeguard)
               lowerB = obj.gamma * (obj.oldEtak^obj.alpha);
               if lowerB > 0.1 && lowerB > etak
                  etak = lowerB;
                  status = 'Safeguard';
               end

               % Final-stage oversolving safeguard (Eisenstat & Walker 1996)
               tauNL = obj.getTauNL();
               etaOversolve = (obj.alphaSafe * tauNL) / max(physNorm, eps);
               if etak < etaOversolve
                  etak = etaOversolve;
                  status = 'Oversolve';
               end

               % Clamp forcing term within [minEtak, maxEtak]
               if etak > obj.maxEtak
                  etak = obj.maxEtak;
                  status = 'CapMax';
               end

               if etak < obj.minEtak
                  targetTol = obj.minEtak;
                  status = 'FloorMin';
               else
                  targetTol = etak;
               end

               % Persist step state for next Newton iteration and retry detection
               obj.oldEtak = targetTol;
               obj.normFxold = physNorm;
               obj.lastNonlinIter = nonlinIter;
               obj.lastPhysNorm = physNorm;
               obj.lastRawEtak = rawEtak;
               obj.lastRatio = ratio;
            else
               % Start from initial forcing term eta0
               ratio = NaN;
               rawEtak = obj.eta0;
               targetTol = obj.eta0;
               status = 'Init';

               obj.oldEtak = obj.eta0;
               obj.normFxold = physNorm;
               obj.lastNonlinIter = nonlinIter;
               obj.lastPhysNorm = physNorm;
               obj.lastRawEtak = rawEtak;
               obj.lastRatio = ratio;
            end

            obj.unscaledTol = targetTol;
            obj.Tol = targetTol;
         end

         % Cache step metadata until post-solve recordStats commits deltaT and relres
         obj.curStep.ratio = ratio;
         obj.curStep.rawEta = rawEtak;
         obj.curStep.eta = targetTol;
         obj.curStep.tolReq = obj.Tol;
         obj.curStep.nonlinIter = nonlinIter;
         obj.curStep.status = status;
      end

      function recordStats(obj, linsolver)
         % Append post-solve Krylov metrics (relres, iters/decade, deltaT) for the last solve
         relres = linsolver.params.lastRelres;
         
         % Measure achieved residual reduction in decades for iters/decade metric
         D = -log10(max(relres, 100 * eps));

         % Actual recording of the statistics
         ns = linsolver.nSolve;
         idx = numel(obj.stats.time) + 1;
         obj.stats.time(idx)       = linsolver.timeLin(ns);
         obj.stats.nonlinIter(idx) = obj.curStep.nonlinIter;
         obj.stats.ratio(idx)      = obj.curStep.ratio;
         obj.stats.rawEta(idx)     = obj.curStep.rawEta;
         obj.stats.eta(idx)        = obj.curStep.eta;
         obj.stats.tolReq(idx)     = obj.curStep.tolReq;
         obj.stats.tolReached(idx) = relres;
         obj.stats.iter(idx)       = linsolver.params.iter;
         obj.stats.itPerDec(idx) = max(1, linsolver.params.iter) / D;
         obj.stats.deltaT(idx)     = linsolver.Delta_T(ns);
         obj.stats.status{idx}     = obj.curStep.status;
      end

      function printStats(obj)
         % Print per-solve Eisenstat-Walker history table and forcing-term summary
         if isempty(obj.stats.time)
            return;
         end

         lineLen = 115;
         sepLine = [repmat('-', 1, lineLen), '\n'];
         titleStr = 'Eisenstat-Walker Convergence Statistics';
         padLeft = floor((lineLen - length(titleStr)) / 2);
         padRight = lineLen - length(titleStr) - padLeft;

         fprintf(['\n', sepLine]);
         fprintf([repmat(' ', 1, padLeft), titleStr, repmat(' ', 1, padRight), '\n']);
         fprintf(sepLine);
         fprintf('| %8s | %5s | %8s | %8s | %8s | %8s | %8s | %5s | %6s | %8s | %-9s |\n', ...
            'PhysTime', 'NL-It', 'Ratio', 'RawEta', 'Eta', 'TolReq', 'TolReach', 'LinIt', 'It/Dec', 'DeltaT', 'Status');
         fprintf(sepLine);

         % Format rows, substituting dashes for NaN metrics (first, linear, or zero-reduction solves)
         n = length(obj.stats.time);
         for i = 1:n
            if isnan(obj.stats.ratio(i))
               ratioStr = '     ---';
            else
               ratioStr = sprintf('%8.2e', obj.stats.ratio(i));
            end

            if isnan(obj.stats.rawEta(i))
               rawEtaStr = '     ---';
            else
               rawEtaStr = sprintf('%8.2e', obj.stats.rawEta(i));
            end

            if isnan(obj.stats.itPerDec(i))
               itPerDecStr = '   ---';
            else
               itPerDecStr = sprintf('%6.1f', obj.stats.itPerDec(i));
            end

            fprintf('| %.2e | %5d | %8s | %8s | %.2e | %.2e | %.2e | %5d | %6s | %.2e | %-9s |\n', ...
               obj.stats.time(i), ...
               obj.stats.nonlinIter(i), ...
               ratioStr, ...
               rawEtaStr, ...
               obj.stats.eta(i), ...
               obj.stats.tolReq(i), ...
               obj.stats.tolReached(i), ...
               obj.stats.iter(i), ...
               itPerDecStr, ...
               obj.stats.deltaT(i), ...
               obj.stats.status{i});
         end
         fprintf(sepLine);
         fprintf('Eisenstat-Walker Summary:\n');
         fprintf('  Average Forcing Term (Eta)    = %e (min: %e, max: %e)\n', ...
            mean(obj.stats.eta), min(obj.stats.eta), max(obj.stats.eta));
         fprintf('  Average Krylov Iters / Decade = %.1f\n', mean(obj.stats.itPerDec, 'omitnan'));

         % Tally occurrences of each forcing-term status
         statusList = obj.stats.status;
         uniqueStatuses = unique(statusList);
         statusCountStr = '';
         for s = 1:length(uniqueStatuses)
            st = uniqueStatuses{s};
            cnt = sum(strcmp(statusList, st));
            if s > 1
               statusCountStr = [statusCountStr, ', ']; %#ok<AGROW>
            end
            statusCountStr = [statusCountStr, sprintf('%s=%d', st, cnt)]; %#ok<AGROW>
         end
         fprintf('  Forcing term status counts    : %s\n', statusCountStr);
         fprintf(sepLine);
         fprintf('\n');
      end
   end
end
