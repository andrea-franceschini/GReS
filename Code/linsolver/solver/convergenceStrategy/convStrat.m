classdef (Abstract) convStrat < handle
   % Abstract class to handle to choose which convergence strategy is to be used.
   
   properties (SetAccess = protected, GetAccess = public)
      % Tolerance value
      Tol

      % Tolerance for linear step
      linearTol
   end

   methods
      function obj = convStrat(generalsolver)
         % Require exactly 1 non-empty argument
         arguments
            generalsolver (1,1) {mustBeNonempty}
         end
         
         % Select the linear tolerance to be the relative tolerance of the
         % nonlinear solver
         obj.linearTol = generalsolver.simparams.relTol;
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