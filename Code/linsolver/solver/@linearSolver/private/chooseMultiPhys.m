function [ChronosFlag,Prec] = chooseMultiPhys(obj,generalsolver,debugflag,physname)

   % List of allowed physics
   allowedPhysics = {'pressure', 'displacements','fractureJump','saturation'};
   
   % Check if any entry in physname is NOT among the allowed physics
   if any(~ismember(physname, allowedPhysics))
       gresLog().warning(3, 'Multiphysics not yet supported');
       if gresLog().getVerbosity() >= 3
           disp([obj.systemInfo.varBlocks.varName]);
       end
       Prec = [];
       ChronosFlag = false;

   % Exactly 2 physics: ('displacements' AND 'pressure')
   elseif numel(unique(physname)) == 2 && ...
          ismember("displacements", physname) && ...
          ismember("pressure", physname)
   
       Prec = fixedStress(debugflag, generalsolver);
       ChronosFlag = true;

   % Exactly 2 physics: ('pressure' AND 'saturation')
   elseif numel(unique(physname)) == 2 && ...
          ismember("pressure", physname) && ...
          ismember("saturation", physname)
   
       Prec = CPR(debugflag, generalsolver);
       ChronosFlag = true;

   % Exactly 2 physics: ('displacements' AND 'fractureJump')
   elseif numel(unique(physname)) == 2 && ...
          ismember("displacements", physname) && ...
          ismember("fractureJump", physname)
   
       Prec = efemPrec(debugflag, generalsolver, obj.nsyTol);
       ChronosFlag = true;

   % Exactly 3 physics: ('displacements' AND 'pressure' AND 'saturation') -> White et al. (2019)
   elseif numel(unique(physname)) == 3 && ...
          ismember("displacements", physname) && ...
          ismember("pressure", physname) && ...
          ismember("saturation", physname)
   
       Prec = fixedStress(debugflag, generalsolver, 'flowPreconditioner', 'CPR');
       ChronosFlag = true;
   
   % Any other combination of allowed physics not explicitly supported
   else
       gresLog().warning(3, 'Multiphysics not yet supported');
       Prec = [];
       ChronosFlag = false;
   end
end
