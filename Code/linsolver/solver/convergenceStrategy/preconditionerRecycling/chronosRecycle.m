classdef chronosRecycle < precRecycle
   % Tolerance-normalized Chronos preconditioner and SAM recycling strategy.
   %
   % Decouples variable linear tolerances from preconditioner degradation
   % tracking by comparing realized solve time against tolerance-scaled
   % baseline work under x0 = 0 and deferred baseline anchoring.

   properties (SetAccess = protected, GetAccess = public)
      % Minimum iteration count required to anchor a baseline solve
      minBaselineIter = 5

      % Residual floor to prevent log10(0) divergence
      minRelresFloor = 100 * eps
   end

   properties (Access = private)
      % Deferred baseline state for preconditioner recycling
      precBaselineReady = false
      precBaselineSkipCount = 0

      % Deferred baseline state for SAM recycling
      samBaselineReady = false
      samBaselineSkipCount = 0
   end

   methods
      function obj = chronosRecycle(generalsolver)
         % Initialize tolerance-normalized Chronos recycling strategy
         obj@precRecycle(generalsolver);
      end

      function recomputePrec(obj, linsolver, Tend)
         % Track tolerance-normalized excess solve time and trigger rebuild when Delta_T > alpha * tSetup
         iter = linsolver.params.iter;
         D = obj.computeResidualReduction(linsolver.params.lastRelres);
         Dref = obj.computeRefReduction(obj.generalsolver.simparams.relTol);

         % Non-positive alpha disables recycling and forces rebuild every solve
         if linsolver.alpha <= 0.
            linsolver.Delta_T(linsolver.nSolve) = 0;
            linsolver.requestPrecComp = true;
            obj.precBaselineReady = false;
            obj.precBaselineSkipCount = 0;
            obj.samBaselineReady = false;
            obj.samBaselineSkipCount = 0;
         elseif linsolver.requestPrecComp
            % Reset counters, clear SAM baseline, and anchor T0 if representative
            linsolver.cumTSolveAfterPrec = 0;
            linsolver.requestPrecComp = false;
            linsolver.Delta_T(linsolver.nSolve) = 0;

            % SAM was reset in SolveLin when the preconditioner was recomputed
            obj.samBaselineReady = false;
            obj.samBaselineSkipCount = 0;

            if iter >= obj.minBaselineIter && D > 0
               linsolver.params.firstSolveTAfterPrecComp = Tend * (Dref / D);
               obj.precBaselineReady = true;
            else
               % Defer baseline anchoring until a representative solve occurs
               obj.precBaselineReady = false;
            end
            obj.precBaselineSkipCount = linsolver.params.nSolveSinceLastPrecComp;
         elseif ~obj.precBaselineReady
            % Wait for a representative solve (iter >= minBaselineIter and D > 0) to anchor T0
            if iter >= obj.minBaselineIter && D > 0
               linsolver.params.firstSolveTAfterPrecComp = Tend * (Dref / D);
               linsolver.cumTSolveAfterPrec = 0;
               obj.precBaselineReady = true;
               obj.precBaselineSkipCount = linsolver.params.nSolveSinceLastPrecComp;
            end
            linsolver.Delta_T(linsolver.nSolve) = 0;
         else
            % Accumulate excess wall seconds vs baseline and test break-even
            T0 = linsolver.params.firstSolveTAfterPrecComp;
            Teff = obj.computeEffectiveTime(Tend, D, Dref, T0, iter);

            linsolver.cumTSolveAfterPrec = linsolver.cumTSolveAfterPrec + Teff;
            nAccum = linsolver.params.nSolveSinceLastPrecComp - obj.precBaselineSkipCount;
            deltaT = linsolver.cumTSolveAfterPrec - nAccum * T0;
            linsolver.Delta_T(linsolver.nSolve) = deltaT;

            % Look up setup cost of the currently active preconditioner
            tSetup = linsolver.precCompLin(end - linsolver.params.nSolveSinceLastPrecComp);
            if deltaT > linsolver.alpha * tSetup
               linsolver.requestPrecComp = true;
            end
         end
      end

      function recomputeSAM(obj, linsolver, sam, Tend)
         % Build initial SAM on mild prec degradation and recycle active SAM via same break-even rule
         if ~linsolver.useSAM
            return
         end

         tSetup = linsolver.precCompLin(end - linsolver.params.nSolveSinceLastPrecComp);
         nsolve = linsolver.nSolve;

         if ~sam.precJustComputed
            % Trigger first SAM build once main prec degrades past firstCompPercDegrad * tSetup
            if linsolver.Delta_T(nsolve) > sam.firstCompPercDegrad * tSetup && isempty(sam.N)
               sam.requestComp = true;
               return
            end

            if ~isempty(sam.N)
               iter = linsolver.params.iter;
               D = obj.computeResidualReduction(linsolver.params.lastRelres);
               Dref = obj.computeRefReduction(obj.generalsolver.simparams.relTol);

               if sam.requestComp
                  % Reset SAM counters and anchor T0SAM if representative
                  sam.cumTSolveAfterComp = 0;
                  sam.Delta_T(nsolve) = 0;
                  sam.requestComp = false;

                  if iter >= obj.minBaselineIter && D > 0
                     sam.firstSolveTAfterComp = Tend * (Dref / D);
                     obj.samBaselineReady = true;
                  else
                     obj.samBaselineReady = false;
                  end
                  obj.samBaselineSkipCount = sam.nSolveSinceLastComp;
               elseif ~obj.samBaselineReady
                  % Deferred SAM anchoring until a representative solve occurs
                  if iter >= obj.minBaselineIter && D > 0
                     sam.firstSolveTAfterComp = Tend * (Dref / D);
                     sam.cumTSolveAfterComp = 0;
                     obj.samBaselineReady = true;
                     obj.samBaselineSkipCount = sam.nSolveSinceLastComp;
                  end
                  sam.Delta_T(nsolve) = 0;
                  if sam.alpha < 0.
                     sam.requestComp = true;
                  end
               else
                  % Accumulate SAM excess time and rebuild when Delta_T > sam.alpha * tSetupSAM
                  T0SAM = sam.firstSolveTAfterComp;
                  TeffSAM = obj.computeEffectiveTime(Tend, D, Dref, T0SAM, iter);

                  sam.cumTSolveAfterComp = sam.cumTSolveAfterComp + TeffSAM;
                  nAccumSAM = sam.nSolveSinceLastComp - obj.samBaselineSkipCount;
                  sam.Delta_T(nsolve) = sam.cumTSolveAfterComp - nAccumSAM * T0SAM;

                  tSetupSAM = sam.CompLin(end - sam.nSolveSinceLastComp);
                  if sam.Delta_T(nsolve) > sam.alpha * tSetupSAM || sam.alpha < 0.
                     sam.requestComp = true;
                  end
               end
            end
         else
            % Main preconditioner was just rebuilt; skip SAM monitoring on this solve
            sam.precJustComputed = false;
            obj.samBaselineReady = false;
            obj.samBaselineSkipCount = 0;
         end
      end
   end

   methods (Access = private)
      function D = computeResidualReduction(obj, relres)
         % Achieved residual reduction in decades: D = -log10(max(relres, floor))
         if isempty(relres) || ~isfinite(relres) || relres >= 1
            D = 0;
            return
         end
         D = -log10(max(relres, obj.minRelresFloor));
      end

      function Dref = computeRefReduction(obj, refTol)
         % Reference reduction in decades corresponding to nonlinear relTol
         tolClamped = min(max(refTol, obj.minRelresFloor), 1 - eps);
         Dref = -log10(tolClamped);
      end

      function Teff = computeEffectiveTime(~, Tend, D, Dref, T0, iter)
         % Excess wall seconds realized over a fresh prec solving D decades (sum(Teff - T0));
         % non-reducing or zero-iteration solves contribute zero excess (Teff = T0).
         if iter <= 0 || D <= 0
            Teff = T0;
         else
            Teff = T0 + (Tend - (D / Dref) * T0);
         end
      end
   end
end
