function [isLinear] = getIsLinear(generalsolver)
% GETISLINEAR Checks if all domain physics solvers and interfaces are linear

   % Initialize to false to return at the first nonlinear solver in any domain
   isLinear = false;

   % Loop over the different domains
   for i = 1:generalsolver.nDom

      % Get current domain handle
      dom = generalsolver.domains(i);

      % Loop over the various solvers of this domain
      solver = dom.solverNames;
      for j = 1:numel(solver)

         % Call isLinear on the current solver in domain i
         lin = dom.getPhysicsSolver(solver(j)).isLinear();
   
         % Early exit, if one domain has a solver which is nonlinear
         % then all the system is nonlinear
         if lin == false
            return;
         end
      end
   end

   % Loop over the different interfaces
   for i = 1:generalsolver.nInterf

      % Call isLinear on the current interface
      lin = generalsolver.interfaces{i}.isLinear();

      % Early exit, if one interface is nonlinear
      % then all the system is nonlinear
      if lin == false
         return
      end      
   end

   % If reached here all the domains solvers are linear
   isLinear = true;
end
