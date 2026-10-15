classdef EisenstatWalker < convStrat
   % Eisenstat-Walker Choice 2 convergence strategy (Eisenstat & Walker, 1996)
   
   properties (SetAccess = private, GetAccess = public)
      % Eisenstat-Walker parameters (Choice 2)
      gamma = 1.0
      alpha = 0.5 * (1 + sqrt(5))
      maxEtak = 0.5
      minEtak
      oldEtak
      eta0 = 0.1
      normFxold
      unscaledTol          % Physical forcing term eta_k (independent of Ruiz scaling)
      
      % Preconditioner recycling tracking
      initNormSolveT = inf % Baseline normalized solve time (sec per decade)
      initConvRate = inf   % Baseline Krylov rate (iters per decade)
      degrad = 1.5
      highPrec = 1e-4
      highDeg = 2.0
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
         obj@convStrat(generalsolver);
         obj.minEtak = generalsolver.simparams.relTol;
         obj.linearTol = generalsolver.simparams.relTol;

         % Parse optional EW parameters from linSolverParams if provided
         if isfield(generalsolver.simparams, 'linSolverParams')
            lsp = generalsolver.simparams.linSolverParams;
            if isfield(lsp, 'eta0')
               obj.eta0 = str2double(string(lsp.eta0));
            end
            if isfield(lsp, 'maxEtak')
               obj.maxEtak = str2double(string(lsp.maxEtak));
            end
            if isfield(lsp, 'gamma')
               obj.gamma = str2double(string(lsp.gamma));
            end
         end
      end

      function computeTol(obj, b, islinear, nonlinIter, Ruiz, backupb)
         % Extract unscaled physical residual norm
         physNorm = norm(cell2matrix(backupb));

         % Check if the solver is retrying to converge
         isRetry = (nonlinIter == obj.lastNonlinIter && abs(physNorm - obj.lastPhysNorm) <= 1e-12 * max(1, physNorm));

         % Compute target physical inexact Newton forcing term eta_k
         if islinear
            targetTol = obj.minEtak;
            rawEtak = obj.minEtak;
            ratio = NaN;
            status = 'Linear';
         elseif isRetry
            % Retry of current Newton step: preserve tolerance without altering history
            targetTol = obj.unscaledTol;
            rawEtak = obj.lastRawEtak;
            ratio = obj.lastRatio;
            status = 'Retry';
         else
            if nonlinIter ~= 1
               ratio = physNorm / obj.normFxold;
               rawEtak = obj.gamma * (ratio)^obj.alpha;
               etak = rawEtak;
               status = 'Choice2';

               % Choice 2 safeguard: prevent etak from dropping too rapidly
               lowerB = obj.gamma * (obj.oldEtak^obj.alpha);
               if lowerB > 0.1 && lowerB > etak
                  etak = lowerB;
                  status = 'Safeguard';
               end

               % Cap at globalization upper bound
               if etak > obj.maxEtak
                  etak = obj.maxEtak;
                  status = 'CapMax';
               end

               % Lower bound from nonlinear solver stopping criteria
               if etak < obj.minEtak
                  targetTol = obj.minEtak;
                  status = 'FloorMin';
               else
                  targetTol = etak;
               end

               % Update state history
               obj.oldEtak = targetTol;
               obj.normFxold = physNorm;
               obj.lastNonlinIter = nonlinIter;
               obj.lastPhysNorm = physNorm;
               obj.lastRawEtak = rawEtak;
               obj.lastRatio = ratio;
            else
               % Initial iteration
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
         end

         obj.unscaledTol = targetTol;
         obj.Tol = targetTol;

         % Adjust tolerance if Ruiz scaling is active
         if Ruiz.scalingFlag
            b_unscaled_norm = physNorm;
            b_scaled_norm = norm(b);
            if b_scaled_norm > 0 && b_unscaled_norm > 0
               minD = min(Ruiz.fullD);
               scaleFactor = (minD * b_unscaled_norm) / b_scaled_norm;
               obj.Tol = obj.Tol * min(1, scaleFactor);
            end
         end

         % Store step info to be committed to history in recomputePrec
         obj.curStep.ratio = ratio;
         obj.curStep.rawEta = rawEtak;
         obj.curStep.eta = targetTol;
         obj.curStep.tolReq = obj.Tol;
         obj.curStep.nonlinIter = nonlinIter;
         obj.curStep.status = status;
      end

      function recomputePrec(obj, linsolver, Tend)
         % Decades of residual reduction requested by Eisenstat-Walker
         decades = max(0.5, -log10(obj.unscaledTol));

         if linsolver.alpha <= 0.
            linsolver.Delta_T(linsolver.nSolve) = 0;
            linsolver.requestPrecComp = true;
            obj.recordStats(linsolver, decades);
            return
         end

         if linsolver.requestPrecComp
            % --- Preconditioner was freshly computed ---
            linsolver.params.firstSolveTAfterPrecComp = Tend;
            linsolver.cumTSolveAfterPrec = 0;
            linsolver.requestPrecComp = false;
            linsolver.Delta_T(linsolver.nSolve) = 0;
        
            % Record normalized baseline metrics per decade
            obj.initNormSolveT = Tend / decades;
            obj.initConvRate = max(1, linsolver.params.iter) / decades;
         else
            % Normalized solve time and Krylov iteration rate for current solve
            currentNormSolveT = Tend / decades;
            currentConvRate = max(1, linsolver.params.iter) / decades;
            
            % Normalized extra time spent due to preconditioner degradation
            normDeltaT = max(0, currentNormSolveT - obj.initNormSolveT) * decades;
            linsolver.cumTSolveAfterPrec = linsolver.cumTSolveAfterPrec + normDeltaT;
            linsolver.Delta_T(linsolver.nSolve) = linsolver.cumTSolveAfterPrec;
           
            tSetup = linsolver.precCompLin(end - linsolver.params.nSolveSinceLastPrecComp);
            convDegradation = currentConvRate / max(1e-12, obj.initConvRate);
           
            % Evaluate recycling conditions
            timeExceeded = (linsolver.Delta_T(linsolver.nSolve) > linsolver.alpha * tSetup);
            rateDegraded = (convDegradation > obj.degrad);
            tightEtaStall = (obj.unscaledTol < obj.highPrec) && (convDegradation > obj.highDeg);
           
            if (timeExceeded && rateDegraded) || tightEtaStall
                linsolver.requestPrecComp = true;
            end
         end

         % Record statistics for this linear solve
         obj.recordStats(linsolver, decades);
      end

      function printStats(obj)
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

            fprintf('| %.2e | %5d | %8s | %8s | %.2e | %.2e | %.2e | %5d | %6.1f | %.2e | %-9s |\n', ...
               obj.stats.time(i), ...
               obj.stats.nonlinIter(i), ...
               ratioStr, ...
               rawEtaStr, ...
               obj.stats.eta(i), ...
               obj.stats.tolReq(i), ...
               obj.stats.tolReached(i), ...
               obj.stats.iter(i), ...
               obj.stats.itPerDec(i), ...
               obj.stats.deltaT(i), ...
               obj.stats.status{i});
         end
         fprintf(sepLine);
         fprintf('Eisenstat-Walker Summary:\n');
         fprintf('  Average Forcing Term (Eta)    = %e (min: %e, max: %e)\n', ...
            mean(obj.stats.eta), min(obj.stats.eta), max(obj.stats.eta));
         fprintf('  Average Krylov Iters / Decade = %.1f\n', mean(obj.stats.itPerDec));
         
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

      function recomputeSAM(obj, linsolver, sam, Tend)
         if linsolver.useSAM
            tSetup = linsolver.precCompLin(end - linsolver.params.nSolveSinceLastPrecComp);
            nsolve = linsolver.nSolve;

            if ~sam.precJustComputed
               if linsolver.Delta_T(nsolve) > sam.firstCompPercDegrad * tSetup && isempty(sam.N)
                  sam.requestComp = true;
                  return;
               end

               if ~isempty(sam.N)
                  if sam.requestComp
                     sam.firstSolveTAfterComp = Tend;
                     sam.cumTSolveAfterComp = 0;
                     sam.Delta_T(nsolve) = 0;
                     sam.requestComp = false;
                  else
                     sam.cumTSolveAfterComp = sam.cumTSolveAfterComp + Tend;
                     sam.Delta_T(nsolve) = sam.cumTSolveAfterComp - sam.nSolveSinceLastComp * sam.firstSolveTAfterComp;
               
                     tSetupSAM = sam.CompLin(end - sam.nSolveSinceLastComp);
                     if sam.Delta_T(nsolve) > sam.alpha * tSetupSAM || sam.alpha < 0.
                        sam.requestComp = true;
                     end
                  end
               end
            else
               sam.precJustComputed = false;
            end
         end
      end
   end

   methods (Access = private)
      function recordStats(obj, linsolver, decades)
         idx = linsolver.nSolve;

         obj.stats.time(idx)       = linsolver.timeLin(idx);
         obj.stats.nonlinIter(idx) = obj.curStep.nonlinIter;
         obj.stats.ratio(idx)      = obj.curStep.ratio;
         obj.stats.rawEta(idx)     = obj.curStep.rawEta;
         obj.stats.eta(idx)        = obj.curStep.eta;
         obj.stats.tolReq(idx)     = obj.curStep.tolReq;
         obj.stats.tolReached(idx) = linsolver.params.lastRelres;
         obj.stats.iter(idx)       = linsolver.params.iter;
         obj.stats.itPerDec(idx)   = max(1, linsolver.params.iter) / decades;
         obj.stats.deltaT(idx)     = linsolver.Delta_T(idx);
         obj.stats.status{idx}     = obj.curStep.status;
      end
   end
end