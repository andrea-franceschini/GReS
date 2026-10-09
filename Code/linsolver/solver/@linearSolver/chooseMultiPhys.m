function [ChronosFlag,Prec] = chooseMultiPhys(obj,generalsolver,debugflag,physname,innerPrec,innerPrecFlux)

   if nargin < 5 || isempty(innerPrec)
      innerPrec = 'amg';
   end
   if nargin < 6 || isempty(innerPrecFlux)
      innerPrecFlux = innerPrec;
   end

   % List of allowed physics
   allowedPhysics = {'pressure', 'displacements','fractureJump'};
   
   % Check if any entry in physname is NOT among the allowed physics
   if any(~ismember(physname, allowedPhysics))
       gresLog().warning(3, 'Multiphysics not yet supported');
       if gresLog().getVerbosity() >= 3
           physNames = arrayfun(@(x) x.dofm.getVariableNames(), generalsolver.domains);
           disp(physNames);
       end
       Prec = [];
       ChronosFlag = false;
   
   % Exactly 2 physics: ('displacements' AND 'pressure')
   elseif numel(unique(physname)) == 2 && ...
          ismember("displacements", physname) && ...
          ismember("pressure", physname)
   
       Prec = fixedStress(debugflag, generalsolver, innerPrec, innerPrecFlux);
       ChronosFlag = true;

   % Exactly 2 physics: ('displacements' AND 'fractureJump')
   elseif numel(unique(physname)) == 2 && ...
          ismember("displacements", physname) && ...
          ismember("fractureJump", physname)
   
       Prec = efemPrec(debugflag, generalsolver, obj.nsyTol, innerPrec);
       ChronosFlag = true;

   % Any other combination of allowed physics not explicitly supported
   else
       gresLog().warning(3, 'Multiphysics not yet supported');
       Prec = [];
       ChronosFlag = false;
   end
end
