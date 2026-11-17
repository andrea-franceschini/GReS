function [Prec,ChronosFlag,scalingFlag] = choosePrec(obj,debugflag,generalsolver,physname)
   
   % Initialize the Chronos and scaling flags
   ChronosFlag = false;
   scalingFlag = false;

   sys = obj.systemInfo;
   multiDomFlag = (generalsolver.nDom > 1);

   % Check if the problem comes from multiphysics (max number of variables in any single domain > 1)
   multiPhysFlag = (max(arrayfun(@(x) x.dofm.getNumberOfVariables(), generalsolver.domains)) > 1);

   % Check if multiphysics multidomain
   if multiPhysFlag && multiDomFlag
      scalingFlag = true;
      order = "racpFirst";
      if isfield(generalsolver.simparams.linSolverParams, 'multiDomPrecOrder')
         order = generalsolver.simparams.linSolverParams.multiDomPrecOrder;
      end
      Prec = resolvePreconditioner(obj, generalsolver, 1:sys.nBlocks, order);
      if ~isempty(Prec)
         ChronosFlag = true;
      else
         gresLog().warning(3, 'Multiphysics with multidomain could not be resolved, falling back to MATLAB solver');
         ChronosFlag = false;
      end
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
      [ChronosFlag,Prec] = chooseMultiPhys(obj,generalsolver,debugflag,physname);
   else
      % Keep only the unique one for the single physics
      physname = unique(physname);

      [ChronosFlag,Prec] = chooseSinglePhys(obj,generalsolver,debugflag,physname);
   end
end
