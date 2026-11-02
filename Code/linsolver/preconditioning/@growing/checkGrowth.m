function x0 = checkGrowth(obj, linsolver, b)
   
   % Set default output  
   x0 = linsolver.x0;

   % Check size change
   if obj.sizeComp ~= 0
      oldProbSize = size(x0,1);
      newProbSize = size(b,1);

      % Check how much it changed
      newDiff = newProbSize - oldProbSize;
      obj.sizeDiff = obj.sizeDiff + newDiff;

      % If changed append zeros at the bottom for the size needed to get to
      % the new size
      if newDiff > 0
         x0 = [x0; zeros(newDiff,1)];
         if newDiff > 0
            gresLog().log(3,'changed size from %d to %d\n',oldProbSize,newProbSize);
         end
      end
   end
end