function [gCoord,gWeight] = gaussValue(nNode,dim)
%Generic function for Gauss point value and Gauss weight value
% dim = element dimension
% nNode = total node number of each element

% Linear hexahedron
if dim == 3
   if nNode == 8
      gCoord = zeros(dim,nNode);
      
      gCoord(1,1) = -0.577350269189626;
      gCoord(1,2) = 0.577350269189626;
      gCoord(1,3) = 0.577350269189626;
      gCoord(1,4) = -0.577350269189626;
      gCoord(1,5) = -0.577350269189626;
      gCoord(1,6) = 0.577350269189626;
      gCoord(1,7) = 0.577350269189626;
      gCoord(1,8) = -0.577350269189626;
      
      gCoord(2,1) = -0.577350269189626;
      gCoord(2,2) = -0.577350269189626;
      gCoord(2,3) = 0.577350269189626;
      gCoord(2,4) = 0.577350269189626;
      gCoord(2,5) = -0.577350269189626;
      gCoord(2,6) = -0.577350269189626;
      gCoord(2,7) = 0.577350269189626;
      gCoord(2,8) = 0.577350269189626;
      
      gCoord(3,2) = -0.577350269189626;
      gCoord(3,1) = -0.577350269189626;
      gCoord(3,3) = -0.577350269189626;
      gCoord(3,4) = -0.577350269189626;
      gCoord(3,5) = 0.577350269189626;
      gCoord(3,6) = 0.577350269189626;
      gCoord(3,7) = 0.577350269189626;
      gCoord(3,8) = 0.577350269189626;
      
      gWeight = ones(nNode,dim);
      
   else
       error('Element not available')
   end
      
   

end

