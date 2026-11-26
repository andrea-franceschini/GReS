function recycleStrat = choosePrecRecycle(generalsolver)
% choosePrecRecycle - Instantiate the preconditioner recycling strategy.

   lsp = generalsolver.simparams.linSolverParams;

   % Default to chronosRecycle when recycleStrategy is omitted
   if ~isfield(lsp, 'recycleStrategy')
      recycleStrat = chronosRecycle(generalsolver);
      return
   end

   % Validate scalar string/char input and dispatch case-insensitively
   val = lsp.recycleStrategy;
   if ~(ischar(val) && isrow(val)) && ~(isstring(val) && isscalar(val))
      error('choosePrecRecycle:InvalidInputType', ...
         'Invalid type or dimension for "recycleStrategy". Expected scalar string or char row vector.');
   end

   stratName = lower(strtrim(string(val)));
   switch stratName
      case "chronos"
         recycleStrat = chronosRecycle(generalsolver);
      case "hypredrive"
         recycleStrat = hypredrive(generalsolver);
      otherwise
         error('choosePrecRecycle:UnknownStrategy', ...
            'Unknown recycleStrategy "%s". Expected "chronos" or "hypredrive".', stratName);
   end
end
