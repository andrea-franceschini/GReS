classdef constantTol < tolStrategy
   % Strictly constant tolerance strategy (GEOS-style) that always passes
   % the prescribed tolerance regardless of islinear, Ruiz scaling, or alphaSafe.

   properties (SetAccess = private, GetAccess = public)
      % Constant relative tolerance passed to the Krylov solver
      tol = 1e-6
   end

   methods
      function obj = constantTol(generalsolver)
         % Initialize base state and read prescribed constant tolerance from linSolverParams.tol
         obj@tolStrategy(generalsolver);

         if isfield(generalsolver.simparams.linSolverParams, 'tol')
            obj.tol = generalsolver.simparams.linSolverParams.tol;
         end
      end

      function computeTol(obj, ~, ~, ~, ~, ~)
         % Pass prescribed tolerance verbatim without linear, Ruiz, or alphaSafe adjustments
         obj.unscaledTol = obj.tol;
         obj.Tol = obj.tol;
      end
   end
end
