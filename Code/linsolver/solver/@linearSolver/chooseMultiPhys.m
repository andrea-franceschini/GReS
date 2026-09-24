function [ChronosFlag,Prec] = chooseMultiPhys(obj,generalsolver,debugflag,physname)

   % List of allowed physics
   allowedPhysics = {'pressure', 'displacements'};
   
   % Check if any entry in physname is NOT among the allowed physics
   if any(~ismember(physname, allowedPhysics))
       gresLog().warning(3, 'Multiphysics not yet supported');
       if gresLog().getVerbosity() >= 3
           physNames = arrayfun(@(x) x.dofm.getVariableNames(), domainin);
           disp(physNames);
       end
       Prec = [];
       ChronosFlag = false;
   
   % Exactly 2 physics: ('displacements' AND 'pressure')
   elseif numel(unique(physname)) == 2 && ...
          ismember("displacements", physname) && ...
          (ismember("pressure", physname))
   
       Prec = fixedStress(debugflag, generalsolver);
       ChronosFlag = true;

   % Any other combination of allowed physics not explicitly supported
   else
       gresLog().warning(3, 'Multiphysics not yet supported');
       Prec = [];
       ChronosFlag = false;
   end
end
