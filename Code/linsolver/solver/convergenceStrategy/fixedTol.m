classdef fixedTol < convStrat
   % Abstract class to handle to choose which convergence strategy is to be used.
   
   properties (SetAccess = private, GetAccess = public)
      % Tolerance for nonlinear step in the case of fixed tolerances
      nonLinTol = 1e-6
   end

   methods

      function obj = fixedTol(generalsolver)
         % Instantiate the superclass
         obj@convStrat(generalsolver);

         if isfield(generalsolver.simparams.linSolverParams, 'tol')
            % If prescribed by the user then set the tolerance
            obj.nonLinTol = generalsolver.simparams.linSolverParams.tol;
         end
      end

      % Function to compute the tolerance for the linear solver
      function computeTol(obj, b, islinear, ~, Ruiz, backupb)
         
         % If the problem is linear then use the tolerance needed from the nonlinear-solver
         if islinear
            obj.Tol = obj.linearTol;

            % Determine effective tolerance for the scaled linear solver to guarantee
            % convergence of the unscaled physical residual
            if Ruiz.scalingFlag
               b_unscaled_norm = norm(cell2matrix(backupb));
               b_scaled_norm = norm(b);
               if b_scaled_norm > 0 && b_unscaled_norm > 0
                  minD = min(Ruiz.fullD);
                  scaleFactor = (minD * b_unscaled_norm) / b_scaled_norm;
                  obj.Tol = obj.Tol * min(1, scaleFactor);
               end
            end
            return;
         else
            % Nonlinear iteration so use default convergence
            obj.Tol = obj.nonLinTol;
         end
      end

      function recomputePrec(obj,linsolver,Tend)
         % If alpha <=0 then always recompute the preconditioner
         if linsolver.alpha <= 0.
            linsolver.Delta_T(linsolver.nSolve) = 0;
            linsolver.requestPrecComp = true;
            return
         end

         % If the preconditioner has just been computed then do not compute it for the next iter
         if linsolver.requestPrecComp
            % --- Preconditioner was just recomputed ---
            linsolver.params.firstSolveTAfterPrecComp = Tend;
            linsolver.cumTSolveAfterPrec = 0;
            linsolver.requestPrecComp = false;
            linsolver.Delta_T(linsolver.nSolve) = 0;
         else
            % --- Preconditioner was not recomputed last iter ---

            % Update stats
            linsolver.cumTSolveAfterPrec = linsolver.cumTSolveAfterPrec + Tend;
            linsolver.Delta_T(linsolver.nSolve) = linsolver.cumTSolveAfterPrec - ...
                linsolver.params.nSolveSinceLastPrecComp * linsolver.params.firstSolveTAfterPrecComp;
           
            tSetup = linsolver.precCompLin(end - linsolver.params.nSolveSinceLastPrecComp);
            
            % Check if the time for the new solves is sufficient to need a
            % preconditioner recomputation
            if linsolver.Delta_T(linsolver.nSolve) > linsolver.alpha*tSetup
               linsolver.requestPrecComp = true;
            end
         end
      end

      function recomputeSAM(obj,linsolver,sam,Tend)
         % Check if SAM needs to be recomputed, do nothing id not asked to
         % use SAM
         if linsolver.useSAM
            % Get the setup time of the preconditioner and the number of
            % solves since the start
            tSetup = linsolver.precCompLin(end-linsolver.params.nSolveSinceLastPrecComp);
            nsolve = linsolver.nSolve;

            % If the preconditioner has been computed in the last iteration
            % supposedly the matrix now should be quite similar so no use
            % in computing the SAM
            if ~sam.precJustComputed
               % The degradation is enought to compute the sam
               if linsolver.Delta_T(nsolve) > sam.firstCompPercDegrad*tSetup && isempty(sam.N)
                  sam.requestComp = true;
                  return;
               end

               if ~isempty(sam.N)
                  if sam.requestComp
                     % Keep in memory the number of iter it did with the correct matrix
                     sam.firstSolveTAfterComp = Tend;
                     sam.cumTSolveAfterComp = 0;
                     sam.Delta_T(nsolve) = 0;
                     sam.requestComp = false; % Reset the request for SAM computation
                  else
                     % Choose if to recompute the SAM
                     sam.cumTSolveAfterComp = sam.cumTSolveAfterComp + Tend;
                     sam.Delta_T(nsolve) = sam.cumTSolveAfterComp - sam.nSolveSinceLastComp*sam.firstSolveTAfterComp;
               
                     tSetupSAM = sam.CompLin(end-sam.nSolveSinceLastComp);
                     if sam.Delta_T(nsolve) > sam.alpha*tSetupSAM || sam.alpha < 0.
                        sam.requestComp = true;
                     end
                  end
               end
            else
               % The preconditioner has alreay been computed a while ago,
               % can check for SAM computation
               sam.precJustComputed = false;
            end
         end
      end
   end
end