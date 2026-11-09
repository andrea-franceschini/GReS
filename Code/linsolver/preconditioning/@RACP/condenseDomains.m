% Condenses the domains and interfaces in a 2x2 matrix
function A = condenseDomains(obj,A)
   nn = size(A,1);
   idxMain = 1:obj.nDom;
   idxSupp = obj.nDom+1:nn;

   % Treat the multiple domains as if they were one and then use RACP
   if numel(A) ~= 4

      A11 = cell2matrix(A(idxMain,idxMain));
      A12 = cell2matrix(A(idxMain,idxSupp));
      A21 = cell2matrix(A(idxSupp,idxMain));
      A22 = cell2matrix(A(idxSupp,idxSupp));

      clear A;

      A = {A11, A12; A21 A22};
   end
end