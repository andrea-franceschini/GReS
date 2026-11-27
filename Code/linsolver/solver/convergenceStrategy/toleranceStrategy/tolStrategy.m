classdef (Abstract) tolStrategy < handle
   % Abstract base class for linear solver tolerance selection strategies.

   properties (SetAccess = protected, GetAccess = public)
      % General solver handle
      generalsolver

      % Effective tolerance passed to the linear solver
      Tol

      % Unscaled physical tolerance before Ruiz adjustment
      unscaledTol

      % Linear-problem relative tolerance inherited from nonlinear solver
      linearTol

      % Safety factor for oversolving protection
      alphaSafe = 0.01
   end

   methods
      function obj = tolStrategy(generalsolver)
         % Bind general solver handle and inherit nonlinear relTol as linear-step default
         arguments
            generalsolver (1,1) {mustBeNonempty}
         end

         obj.generalsolver = generalsolver;
         obj.linearTol = generalsolver.simparams.relTol;
      end

      function tauNL = getTauNL(obj)
         % Compute effective nonlinear stopping target max(absTol, relTol * ||r_0||)
         tauNL = obj.generalsolver.simparams.absTol;
         tauNL = max(tauNL, obj.generalsolver.simparams.relTol * obj.generalsolver.rhsNormIt0);
      end

      function recordStats(~, ~)
         % Default no-op post-solve hook for strategies without per-solve logs
      end

      function printStats(~)
         % Default no-op print hook for strategies without per-solve tables
      end
   end

   methods (Access = protected)
      function physNorm = getPhysNorm(~, backupb)
         % Return Euclidean norm of the unscaled physical RHS
         physNorm = norm(cell2matrix(backupb));
      end

      function applyLinearRuizTol(obj, targetTol, b, Ruiz, backupb, physNorm)
         % Set target tolerance and tighten scaled tolerance when Ruiz
         % equilibration is active to preserve unscaled residual convergence.
         obj.unscaledTol = targetTol;
         obj.Tol = targetTol;

         if Ruiz.scalingFlag
            if nargin < 6 || isempty(physNorm)
               b_unscaled_norm = obj.getPhysNorm(backupb);
            else
               b_unscaled_norm = physNorm;
            end
            b_scaled_norm = norm(b);
            if b_scaled_norm > 0 && b_unscaled_norm > 0
               % Tighten Tol via min(D) so scaled convergence implies ||r_phys|| <= targetTol * ||b_phys||
               minD = min(Ruiz.fullD);
               scaleFactor = (minD * b_unscaled_norm) / b_scaled_norm;
               obj.Tol = max(obj.Tol * min(1, scaleFactor), 100 * eps);
            end
         end
      end
   end

   methods (Abstract)
      % Populate obj.unscaledTol and obj.Tol before the Krylov solve
      computeTol(obj, b, islinear, nonlinIter, Ruiz, backupb);
   end
end
