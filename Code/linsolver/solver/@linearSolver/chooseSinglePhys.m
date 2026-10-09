function [ChronosFlag,Prec] = chooseSinglePhys(obj,generalsolver,debugflag,physname,precType,innerPrec)

   if nargin < 5 || isempty(precType)
      precType = 'amg';
   end
   if nargin < 6 || isempty(innerPrec)
      innerPrec = 'amg';
   end

   nInt = generalsolver.nInterf;

   % List of allowed physics
   allowedPhysics = {'pressure', 'u', 'displacements'};

   % Supported Single Physics
   if contains(physname, allowedPhysics)
      
      if nInt == 0
         % No interface, its a simple single domain single physics problem.
         % Choose preconditioner via precType flag ('amg' or 'fsai', defaults to 'amg')
         Prec = createInnerPrec(precType,debugflag,generalsolver,physname);
      else
         % Multiple interfaces mean multiple blocks, RACP is needed with chosen inner preconditioner
         Prec = RACP(debugflag,generalsolver,physname,innerPrec);
      end

      ChronosFlag = true;
   else
      % Not a supported physics
      if debugflag
         warning('No preconditioner available for this physics, falling back to matlab solver');
         disp(physname);
      end

      Prec = [];
      ChronosFlag = false;
   end
end
