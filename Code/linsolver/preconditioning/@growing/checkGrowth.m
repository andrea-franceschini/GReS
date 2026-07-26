function x0 = checkGrowth(obj, linsolver, b)
   
   % Set default output  
   x0 = linsolver.x0;

   % Check size change
   oldProbSize = size(x0,1);
   if oldProbSize ~= 0
      newProbSize = size(b,1);

      % Check how much it changed
      obj.sizeDiff = obj.sizeDiff + newProbSize - oldProbSize;

      % If changed append zeros at the bottom for the size needed to get to
      % the new size
      if obj.sizeDiff > 0
         x0 = [x0; zeros(obj.sizeDiff,1)];
         if newProbSize - oldProbSize > 0
            gresLog().log(3,'changed size from %d to %d\n',oldProbSize,newProbSize);
         end
      end
   end
end