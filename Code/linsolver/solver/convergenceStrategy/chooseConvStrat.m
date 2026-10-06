
% Function to choose which convergence strategy is to be used in the
% simulation
function [convergenceStrat] = chooseConvStrat(generalsolver)

   useEW = false;

   % Check if Eisenstat-Walker option is enabled in params
   if isfield(generalsolver.simparams.linSolverParams, 'useEW')
      useEW = generalsolver.simparams.linSolverParams.useEW;
   end

   if useEW
      % Use Eisenstat-Walker
      convergenceStrat = EisenstatWalker(generalsolver);
   else
      % Use fixed nonlinear tolerance
      convergenceStrat = fixedTol(generalsolver);
   end
end