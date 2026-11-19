% Function to choose which convergence strategy is to be used in the simulation
function [convergenceStrat] = chooseConvStrat(generalsolver)
   useEW = false;
   
   % Check if Eisenstat-Walker option is specified in params
   if isfield(generalsolver.simparams.linSolverParams, 'useEW')
      val = generalsolver.simparams.linSolverParams.useEW;
   
      if isnumeric(val) && isscalar(val) && (val == 0 || val == 1)
         useEW = logical(val);
      else
         error('chooseConvStrat:InvalidInputType', ...
            'Invalid type or dimension for "useEW". Expected 0/1 numeric.');
      end
   end
   
   if useEW
      % Use Eisenstat-Walker
      convergenceStrat = EisenstatWalker(generalsolver);
   else
      % Use fixed nonlinear tolerance
      convergenceStrat = fixedTol(generalsolver);
   end
end