function y = applyMultiRACP(obj,A11_aug,B1,B2,inv_D22,x)
   
   % Get the size of the System from the rhs
   nn = size(x,1);

   % Get the size of the inteface blocks from the last block
   sizeInt = size(inv_D22,1);

   % Start size
   n0 = 0;

   % Initialize
   y = x;

   % Get the last block rhs
   x2 = x(nn-sizeInt+1:end,:);

   % Compute the first part
   b1 = x(1:nn-sizeInt) + B1*inv_D22*x2;
   
   for i = 1:obj.nDom
      % Get the size of this domain
      sizeDomI = size(A11_aug{i},1);

      % Apply the AMG of the current domain
      y(n0+1:n0+sizeDomI,1) = obj.AMG{i}.Apply_L(b1(n0+1:n0+sizeDomI));

      % Update the starting size for x
      n0 = n0 + sizeDomI;
   end

   % Complete the RACP application of the current domain
   y(nn-sizeInt+1:end) = inv_D22*(B2*y(1:nn-sizeInt,1) - x2);

end