function [Prec,ChronosFlag] = choosePrec(obj,debugflag,generalsolver,physname)
   
   % Initialize the Chronos Flag
   ChronosFlag = false;

   domainin = generalsolver.domains;

   % Check if the problem comes from multiphysics
   multiPhysFlag = false;
   if(domainin(1).dofm.getNumberOfVariables() > 1) && isempty(physname)
      multiPhysFlag = true; 
      Prec = [];
      gresLog().warning(3,'Multiphysics not yet supported');
      if gresLog().getVerbosity() >= 3
         disp(domainin(1).dofm.getVariableNames());
      end
      % return
   end

   nInt = generalsolver.nInterf;

   % Select the physics, check if asked by user directly
   if isempty(physname)
      physname = domainin(1).dofm.getVariableNames();
   end

   % Needs the growing preconditioner
   if contains(domainin(1).solverNames,"Sedimentation") && ~multiPhysFlag
      Prec = growing(debugflag,generalsolver,physname);
      ChronosFlag = true;
      obj.useSAM = false;
      obj.SAM = [];
      return;
   end

   if ~multiPhysFlag
      % Supported Single Physics
      if contains(physname, {'pressure', 'u','displacements'})
         
         if nInt == 0
            % No interface, its a simple single domain single physics problem,
            % AMG handles it beautifully
            Prec = aAMG(debugflag,generalsolver,physname);
         else
            % Multiple interfaces mean multiple blocks, RACP is needed
            Prec = RACP(debugflag,generalsolver,physname);
         end
      else
         if debugflag
            warning('No preconditioner available for this physics, falling back to matlab solver');
            disp(physname);
         end
         return
      end
      ChronosFlag = true;
   else
      % physname
      if any(contains(physname, {'pressure', 'u'})) && any(contains(physname, 'displacements'))
         Prec = fixedStress(debugflag,generalsolver);
         ChronosFlag = true;
      end
   end
end