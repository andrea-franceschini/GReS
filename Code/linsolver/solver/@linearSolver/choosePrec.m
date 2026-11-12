function [Prec,ChronosFlag,scalingFlag] = choosePrec(obj,debugflag,generalsolver,physname)
   
   % Initialize the Chronos and scaling flags
   ChronosFlag = false;
   scalingFlag = false;

   sys = obj.systemInfo;
   multiDomFlag = (generalsolver.nDom > 1);

   % Check if the problem comes from multiphysics (max number of variables in any single domain > 1)
   multiPhysFlag = (max(arrayfun(@(x) x.dofm.getNumberOfVariables(), generalsolver.domains)) > 1);

   % Early exit with multiphysics multidomain
   if multiPhysFlag && multiDomFlag
      scalingFlag = true;
      gresLog().warning(3,'Multiphysics with multidomain not yet supported');
      Prec = [];
      return;
   end

   % Select the physics, check if asked by user directly
   if isempty(physname)
      physname = [sys.varBlocks.varName];
   else
      multiPhysFlag = numel(physname) > 1;
   end

   % Check if it needs the growing preconditioner
   solvers = [sys.varBlocks.solvers];
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
      scalingFlag = true;
      [ChronosFlag,Prec] = obj.chooseMultiPhys(generalsolver,debugflag,physname);
   else
      % Keep only the unique one for the single physics
      physname = unique(physname);

      [ChronosFlag,Prec] = obj.chooseSinglePhys(generalsolver,debugflag,physname);
   end
end
