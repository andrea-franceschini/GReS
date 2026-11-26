function tolStrat = chooseTolStrategy(generalsolver)
% chooseTolStrategy - Instantiate the linear solver tolerance strategy.

   lsp = generalsolver.simparams.linSolverParams;

   % Default to fixedTol when toleranceStrategy is omitted
   if ~isfield(lsp, 'toleranceStrategy')
      tolStrat = fixedTol(generalsolver);
      return
   end

   % Validate scalar string/char input and dispatch case-insensitively
   val = lsp.toleranceStrategy;
   if ~(ischar(val) && isrow(val)) && ~(isstring(val) && isscalar(val))
      error('chooseTolStrategy:InvalidInputType', ...
         'Invalid type or dimension for "toleranceStrategy". Expected scalar string or char row vector.');
   end

   stratName = lower(strtrim(string(val)));
   switch stratName
      case "fixedtol"
         tolStrat = fixedTol(generalsolver);
      case "constanttol"
         tolStrat = constantTol(generalsolver);
      case {"eisenstatwalker", "ew"}
         tolStrat = EisenstatWalker(generalsolver);
      otherwise
         error('chooseTolStrategy:UnknownStrategy', ...
            'Unknown toleranceStrategy "%s". Expected "fixedTol", "constantTol", or "EisenstatWalker".', stratName);
   end
end
