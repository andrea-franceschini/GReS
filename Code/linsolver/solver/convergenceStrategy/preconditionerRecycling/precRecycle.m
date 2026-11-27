classdef (Abstract) precRecycle < handle
   % Abstract base class for preconditioner and SAM recycling strategies.

   properties (SetAccess = protected, GetAccess = public)
      % General solver handle
      generalsolver
   end

   methods
      function obj = precRecycle(generalsolver)
         % Store general solver handle for accessing simulation tolerances and settings
         arguments
            generalsolver (1,1) {mustBeNonempty}
         end

         obj.generalsolver = generalsolver;
      end

      function printStats(~)
         % Default no-op print hook for strategies without per-solve tables
      end
   end

   methods (Abstract)
      % Decide after each solve whether to rebuild the preconditioner or SAM next step
      recomputePrec(obj, linsolver, Tend);
      recomputeSAM(obj, linsolver, sam, Tend);
   end
end
