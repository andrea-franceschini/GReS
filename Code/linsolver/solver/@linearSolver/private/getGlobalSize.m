function [gSize] = getGlobalSize(A)
% GETGLOBALSIZE Calculates total dimension across diagonal/off-diagonal blocks in cell matrix A

   [nBlockRows, nBlockCols] = size(A);
   rowSizes = zeros(nBlockRows, 1);
   
   for i = 1:nBlockRows
      if ~isempty(A{i,i})
         rowSizes(i) = size(A{i,i}, 1);
      else
         for j = 1:nBlockCols
            C = A{i,j};
            if ~isempty(C)
               rowSizes(i) = size(C, 1);
               break;
            end
         end
      end
   end

   gSize = sum(rowSizes);
end
