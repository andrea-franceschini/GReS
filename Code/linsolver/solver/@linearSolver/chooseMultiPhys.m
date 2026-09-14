function [ChronosFlag,Prec] = chooseMultiPhys(obj,generalsolver,debugflag,physname)

   % If the physics is not with pressure, u or displacements is not yet
   % supported
   if any(~contains(physname, {'pressure', 'u', 'displacements'}))
      gresLog().warning(3,'Multiphysics not yet supported');
      if gresLog().getVerbosity() >= 3
         physNames = arrayfun(@(x) x.dofm.getVariableNames(), domainin);
         disp(physNames);
      end
   end

   if any(contains(physname, {'pressure', 'u'})) && any(contains(physname, 'displacements'))
      Prec = fixedStress(debugflag,generalsolver);
      ChronosFlag = true;
   else
      Prec = [];
      ChronosFlag = false;
   end
end