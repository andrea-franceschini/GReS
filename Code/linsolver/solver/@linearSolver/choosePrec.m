function [Prec,ChronosFlag,scalingFlag] = choosePrec(obj,debugflag,generalsolver,physname,precTypeIn,innerPrecIn,innerPrecFluxIn)
   
   % Initialize the Chronos and scaling flags
   ChronosFlag = false;
   scalingFlag = false;

   % Defaults
   precType      = 'amg';
   innerPrec     = 'amg';
   innerPrecFlux = [];

   % Read from generalsolver.simparams.linSolverParams if available
   lsp = generalsolver.simparams.linSolverParams;
   if isfield(lsp, 'innerPrec'),     innerPrec     = lsp.innerPrec;     end
   if isfield(lsp, 'innerPrecFlux'), innerPrecFlux = lsp.innerPrecFlux; end

   % Explicit argument overrides
   if nargin >= 5 && ~isempty(precTypeIn),      precType      = precTypeIn;      end
   if nargin >= 6 && ~isempty(innerPrecIn),     innerPrec     = innerPrecIn;     end
   if nargin >= 7 && ~isempty(innerPrecFluxIn), innerPrecFlux = innerPrecFluxIn; end

   % If innerPrecFlux is not specified, default to innerPrec
   if isempty(innerPrecFlux)
      innerPrecFlux = innerPrec;
   end

   % Get the domains
   domainin = generalsolver.domains;
   multiDomFlag = (generalsolver.nDom > 1);

   % Check if the problem comes from multiphysics
   multiPhysFlag = (max(arrayfun(@(x) x.dofm.getNumberOfVariables(), domainin)) > 1);

   % Early exit with multiphysics multidomain
   if multiPhysFlag && multiDomFlag
      scalingFlag = true;
      gresLog().warning(3,'Multiphysics with multidomain not yet supported');
      Prec = [];
      return;
   end

   % Select the physics, check if asked by user directly
   if isempty(physname)
      physname = arrayfun(@(x) x.dofm.getVariableNames(), domainin, 'UniformOutput', false);
      physname = [physname{:}];
   else
      multiPhysFlag = numel(physname) > 1;
   end
   physname = string(physname);

   % Check if it needs the growing preconditioner
   solvers = unique(arrayfun(@(x) x.solverNames, domainin));
   if contains(solvers,"Sedimentation")
      if ~multiPhysFlag && ~multiDomFlag
         Prec = growing(debugflag,generalsolver,physname,innerPrec);
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
      [ChronosFlag,Prec] = obj.chooseMultiPhys(generalsolver,debugflag,physname,innerPrec,innerPrecFlux);
   else
      % Keep only the unique one for the single physics
      physname = unique(physname);

      [ChronosFlag,Prec] = obj.chooseSinglePhys(generalsolver,debugflag,physname,precType,innerPrec);
   end
end
