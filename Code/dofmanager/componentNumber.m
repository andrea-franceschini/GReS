function out = componentNumber(mesh,in)
   switch in
      case {'Poromechanics'}
         out = mesh.nDim;
      otherwise
         out = 1;
   end
end

