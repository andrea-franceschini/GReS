function [ChronosFlag,Prec] = chooseSinglePhys(obj,generalsolver,debugflag,physname)

   % List of allowed physics
   allowedPhysics = {'pressure', 'u', 'displacements'};

   % Supported Single Physics
   if contains(physname, allowedPhysics)
      
      sys = obj.systemInfo;
      if sys.nInterfaces == 0
         % No interface, its a simple single domain single physics problem,
         % AMG handles it beautifully
         Prec = aAMG(debugflag,generalsolver,physname);
      else
         if all(contains([sys.interfBlocks.interfType], 'MeshTyingTPFA'))
            Prec = multiDomAMG(debugflag,generalsolver,physname);
         else
            % Multiple interfaces mean multiple blocks, RACP is needed
            Prec = RACP(debugflag,generalsolver,physname);
         end
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
