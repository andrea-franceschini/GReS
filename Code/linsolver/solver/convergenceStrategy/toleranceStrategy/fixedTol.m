classdef fixedTol < tolStrategy
   % Fixed nonlinear tolerance strategy with linear Ruiz scaling and
   % optional alphaSafe oversolving protection.

   properties (SetAccess = private, GetAccess = public)
      % Tolerance for nonlinear step
      nonLinTol = 1e-6
   end

   methods
      function obj = fixedTol(generalsolver)
         % Initialize base state and override nonLinTol from linSolverParams.tol if given
         obj@tolStrategy(generalsolver);

         if isfield(generalsolver.simparams.linSolverParams, 'tol')
            obj.nonLinTol = generalsolver.simparams.linSolverParams.tol;
         end
      end

      function computeTol(obj, b, islinear, ~, Ruiz, backupb)
         % Use Ruiz-adjusted linearTol on linear steps; apply alphaSafe relaxation on nonlinear steps
         if islinear
            % Solve to nonlinear relTol with Ruiz scaling protection
            obj.applyLinearRuizTol(obj.linearTol, b, Ruiz, backupb);
         else
            physNorm = obj.getPhysNorm(backupb);
            tauNL = obj.getTauNL();
            tolFromAbs = (obj.alphaSafe * tauNL) / max(physNorm, eps);

            % Relax tolerance if physNorm is already small, capped at maxTol
            maxTol = 0.1;
            targetTol = min(maxTol, max(obj.nonLinTol, tolFromAbs));

            obj.unscaledTol = targetTol;
            obj.Tol = targetTol;
         end
      end
   end
end
