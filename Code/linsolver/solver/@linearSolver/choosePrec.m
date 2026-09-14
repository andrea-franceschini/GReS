function [Prec,ChronosFlag] = choosePrec(obj,debugflag,generalsolver,physname)
   
   % Initialize the Chronos Flag
   ChronosFlag = false;

   % Get the domains
   domainin = generalsolver.domains;
   multiDomFlag = (generalsolver.nDom > 1);

   % Check if the problem comes from multiphysics
   multiPhysFlag = (max(arrayfun(@(x) x.dofm.getNumberOfVariables(), domainin)) > 1);

   % Early exit with multiphysics multidomain
   if multiPhysFlag && multiDomFlag
      gresLog().warning(3,'Multiphysics with multidomain not yet supported');
      Prec = [];
      return;
   end

   % Select the physics, check if asked by user directly
   if isempty(physname)
      physname = arrayfun(@(x) x.dofm.getVariableNames(), domainin, 'UniformOutput', false);
      physname = [physname{:}];
   end

   % Check if it needs the growing preconditioner
   solvers = unique(arrayfun(@(x) x.solverNames, domainin));
   if contains(solvers,"Sedimentation")
      if ~multiPhysFlag && ~multiDomFlag
         Prec = growing(debugflag,generalsolver,physname);
         ChronosFlag = true;
         obj.useSAM = false;
         obj.SAM = [];
         return;
      else
         gresLog().warning(3,'Multiphysics growth and Multidomain growth not yet supported');
         Prec = [];
         return;
      end         
   end

   % Now choose the correct preconditioner for the correct case
   if(multiPhysFlag && ~multiDomFlag)
      [ChronosFlag,Prec] = obj.chooseMultiPhys(generalsolver,debugflag,physname);
   else
      % Keep only the unique one for the single physics
      physname = unique(physname);

      [ChronosFlag,Prec] = obj.chooseSinglePhys(generalsolver,debugflag,physname);
   end
end