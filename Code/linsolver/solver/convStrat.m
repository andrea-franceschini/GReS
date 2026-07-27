classdef convStrat < handle
   % Class to handle switching between constant tolerance and
   % Eisenstat-Walker tolerance in the solution for the linear system

   properties (SetAccess = private, GetAccess = public)
      % Set default to constant value tolerance
      useEW = false

      % Tolerance value
      Tol

      % Parameters for Eisenstat-Walker
      gamma = 1
      alpha = 0.5 * (1 + sqrt(5))
      maxEtak = 0.9
      minEtak
      oldEtak
      eta0 = 0.1
      normFxold
      initConvRate = inf
      degrad = 1.5
      highPrec = 1e-2
      highDeg = 2.0

   end

   methods

      function obj = convStrat(generalsolver)
         % Check if Eisenstat-Walker option is enabled in params
         if isfield(generalsolver.simparams.linSolverParams, 'useEW')
            obj.useEW = generalsolver.simparams.linSolverParams.useEW;
         end

         % Select the minimum Etak to be the relative tolerance of the
         % nonlinear solver
         obj.minEtak = generalsolver.simparams.relTol;

         if ~obj.useEW
            if isfield(generalsolver.simparams.linSolverParams, 'tol')
               % If prescribed by the user then set the tolerance
               obj.Tol = generalsolver.simparams.linSolverParams.tol;
            else
               % Use relative tolerance of the nonlinear solver as default
               % tolerance 
               obj.Tol = generalsolver.simparams.relTol;
            end
         end
      end

      % Function to compute the tolerance for the linear solver
      function computeTol(obj,nonlinIter,b,islinear)
         % If the problem is linear then use the tolerance needed from the nonlinearsolver
         if islinear
            obj.Tol = obj.minEtak;
            return;
         end

         if obj.useEW
            if nonlinIter ~= 1
               newNorm = norm(b);
               etak = obj.gamma*(newNorm/obj.normFxold)^obj.alpha;
   
               % Check upperbound
               etak = min(etak,obj.maxEtak);
   
               lowerB = obj.gamma*(obj.oldEtak^obj.alpha);
               % Check lowerbound
               if lowerB > 0.1
                  obj.Tol = max(etak,lowerB);
               else
                  obj.Tol = etak;
               end

               % Check it does not go under the minimum tolerance needed
               % from the nonlinear solver
               obj.Tol = max(obj.Tol, obj.minEtak);

               % Update old norm
               obj.normFxold = newNorm;
            else
               % First value at the entrance
               obj.Tol = obj.eta0;
               obj.oldEtak = obj.eta0;
               obj.normFxold = norm(b);
            end
         end
      end

      function recomputePrec(obj,linsolver,Tend)
         % If the preconditioner has just been computed then do not compute it for the next iter
         if linsolver.requestPrecComp
            % --- Preconditioner was just recomputed ---
            linsolver.params.firstSolveTAfterPrecComp = Tend;
            linsolver.cumTSolveAfterPrec = 0;
            linsolver.requestPrecComp = false;
            linsolver.Delta_T(linsolver.nSolve) = 0;
        
            if obj.useEW
               % Track baseline Krylov efficiency right after fresh preconditioner setup
               % eta_k is the current Eisenstat-Walker forcing term for this step
               decadeReduction = max(1e-12, -log10(obj.Tol));
               obj.initConvRate = linsolver.params.iter / decadeReduction;
            end
         
         else
            linsolver.cumTSolveAfterPrec = linsolver.cumTSolveAfterPrec + Tend;
            linsolver.Delta_T(linsolver.nSolve) = linsolver.cumTSolveAfterPrec - ...
                linsolver.params.nSolveSinceLastPrecComp * linsolver.params.firstSolveTAfterPrecComp;
           
            tSetup = linsolver.precCompLin(end - linsolver.params.nSolveSinceLastPrecComp);
            
            % --- Evaluate whether to recompute preconditioner ---
            if ~obj.useEW
               % For simple fixed convergence tolerance check if the time
               % for the new solves is sufficient to need a preconditioner
               % recomputation
               if linsolver.Delta_T(linsolver.nSolve) > linsolver.alpha*tSetup || linsolver.alpha < 0.
                  linsolver.requestPrecComp = true;
               end
            else
               
               % Measure current Krylov solver efficiency under EW forcing term
               decadeReduction = max(1e-12, -log10(obj.Tol));
               currentConvRate = linsolver.params.iter / decadeReduction;
            
               % Efficiency degradation factor compared to fresh preconditioner baseline
               convDegradation = currentConvRate / obj.initConvRate;
            
               % Standard time budget exceeded AND Krylov rate degraded (>1.5x slower per decade)
               timeExceeded = (linsolver.Delta_T(linsolver.nSolve) > linsolver.alpha * tSetup) || (linsolver.alpha < 0);
               rateDegraded = (convDegradation > obj.degrad);
            
               % Emergency trigger if EW demands high precision 
               % and the preconditioner is stalling
               tightEtaStall = (obj.Tol < obj.highPrec) && (convDegradation > obj.highDeg);
            
               if (timeExceeded && rateDegraded) || tightEtaStall
                   linsolver.requestPrecComp = true;
               end
            end
         end
      end
   end
end