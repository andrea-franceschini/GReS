classdef (Abstract) convStrat < handle
   % Abstract class to handle to choose which convergence strategy is to be used.
   
   properties (SetAccess = protected, GetAccess = public)
      % General solver handle
      generalsolver

      % Tolerance value
      Tol

      % Tolerance for linear step
      linearTol

      % Safety factor for oversolving protection
      alphaSafe = 0.1
   end

   methods
      function obj = convStrat(generalsolver)
         % Require exactly 1 non-empty argument
         arguments
            generalsolver (1,1) {mustBeNonempty}
         end
         
         obj.generalsolver = generalsolver;

         % Select the linear tolerance to be the relative tolerance of the
         % nonlinear solver
         obj.linearTol = generalsolver.simparams.relTol;
      end

      function tauNL = getTauNL(obj)
         tauNL = obj.generalsolver.simparams.absTol;
         if isprop(obj.generalsolver, 'rhsNormIt0') && ~isempty(obj.generalsolver.rhsNormIt0) && obj.generalsolver.rhsNormIt0 > 0
            tauNL = max(tauNL, obj.generalsolver.simparams.relTol * obj.generalsolver.rhsNormIt0);
         end
      end

      function printStats(obj)
         % Default empty method for convergence strategies that do not
         % track specialized statistics
      end
   end

   methods (Abstract)
      % Function to compute the tolerance for the linear solver
      computeTol(obj, b, islinear, nonlinIter, Ruiz, backupb);

      % Function to check if the preconditioner needs to be recomputed
      recomputePrec(obj,linsolver,Tend);

      % Function to check if the SAM needs to be recomputed
      recomputeSAM(obj,linsolver,sam,Tend);
   end
end