function [ChronosFlag,Prec] = chooseSinglePhys(obj,generalsolver,debugflag,physname,activeBlocks)

   if nargin < 5 || isempty(activeBlocks)
      scopedBlocks = obj.systemInfo.blocks;
   else
      scopedBlocks = obj.systemInfo.blocks(activeBlocks);
   end

   % List of allowed physics
   allowedPhysics = {'pressure', 'u', 'displacements'};

   % Supported Single Physics
   if contains(physname, allowedPhysics)
      
      interfScoped = scopedBlocks([scopedBlocks.kind] == "interface");
      if isempty(interfScoped)
         Prec = aAMG(debugflag,generalsolver,physname);
      else
         % Saddle-point interfaces present, RACP is needed
         Prec = RACP(debugflag,generalsolver,physname);
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
