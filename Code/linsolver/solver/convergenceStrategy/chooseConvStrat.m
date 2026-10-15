% Function to choose which convergence strategy is to be used in the simulation
function [convergenceStrat] = chooseConvStrat(generalsolver)
   useEW = false;
   
   % Check if Eisenstat-Walker option is specified in params
   if isfield(generalsolver.simparams.linSolverParams, 'useEW')
      val = generalsolver.simparams.linSolverParams.useEW;
   
      if islogical(val) && isscalar(val)
         useEW = val;
      elseif isnumeric(val) && isscalar(val) && (val == 0 || val == 1)
         useEW = logical(val);
      elseif (ischar(val) || isstring(val)) && isscalar(val)
         cleanVal = strtrim(string(val));
         if strcmpi(cleanVal, "true") || cleanVal == "1"
            useEW = true;
         elseif strcmpi(cleanVal, "false") || cleanVal == "0"
            useEW = false;
         else
            error('chooseConvStrat:InvalidInput', ...
               'Expected "useEW" to be true/false or 1/0, but received: %s', cleanVal);
         end
      else
         error('chooseConvStrat:InvalidInputType', ...
            'Invalid type or dimension for "useEW". Expected scalar boolean, 0/1 numeric, or valid string/char.');
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