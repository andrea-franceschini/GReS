% Condenses the domains and interfaces in a 2x2 matrix
function A = condenseDomains(obj,A)
   
   % Checks if any of the matrix entries are 0x0 blocks
   [ZeroSpRow,ZeroSpCol] = find(cellfun(@(x) isempty(x), A));
   if ~isempty(ZeroSpRow)
      % get the correct number of rows
      rows = zeros(size(A,1),1);
      for i = 1:size(A,1)
         % Loop over the 
         for j = 1:size(A,1)
            sizz = size(A{i,j},1);
            if sizz ~= 0
               rows(i) = sizz;
               break;
            end
            if j == size(A,1)
               error("no full blocks in this matrix at one column");
            end
         end
      end 
      % assign the correct dimension to the matrices
      for i = 1:length(ZeroSpRow)
         A{ZeroSpRow(i),ZeroSpCol(i)} = sparse(rows(ZeroSpRow(i)),rows(ZeroSpCol(i)));
      end
   end

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