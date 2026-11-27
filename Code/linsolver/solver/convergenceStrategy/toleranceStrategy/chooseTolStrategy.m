function tolStrat = chooseTolStrategy(generalsolver)
% chooseTolStrategy - Instantiate the linear solver tolerance strategy.

   lsp = generalsolver.simparams.linSolverParams;

   if isfield(lsp, 'useEW')
      val = lsp.useEW;
      if isnumeric(val) && isscalar(val) && val == 1
         tolStrat = EisenstatWalker(generalsolver);
         return
      elseif ~(isnumeric(val) && isscalar(val) && val == 0)
         warning('chooseTolStrategy:InvalidInput', ...
            'Expected "useEW" to be 0 or 1; defaulting to fixedTol.');
      end
   end

   tolStrat = fixedTol(generalsolver);
end
