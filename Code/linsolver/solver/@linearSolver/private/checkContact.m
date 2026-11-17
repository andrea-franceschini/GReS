function checkContact(obj)
   if isequal(obj.Prec.phys, 1.1)
      if obj.generalsolver.iterConfig > obj.iterConfigOld && obj.iterConfigOld == 1
         % If the SAM is being used then it is better to just compute it
         % wrt computing the whole preconditioner
         if ~obj.useSAM
            obj.requestPrecComp = true;
         else
            obj.SAM.requestComp = true;
         end
         obj.iterConfigOld = obj.generalsolver.iterConfig;
      elseif obj.generalsolver.iterConfig < obj.iterConfigOld
         obj.iterConfigOld = obj.generalsolver.iterConfig;
      end
   end
end
