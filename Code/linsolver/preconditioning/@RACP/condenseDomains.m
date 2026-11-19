% Condenses only the saddle-point interfaces into block K+1 while keeping variable blocks 1..K separate
function A_out = condenseDomains(obj,A)
   idxMain = obj.getKindBlocks("variable");
   idxSupp = obj.getKindBlocks("interface");
   K = numel(idxMain);
   
   % If already condensed to (K+1) x (K+1) with a single merged interface block
   if size(A, 1) == (K + 1) && size(A, 2) == (K + 1)
      A_out = A;
      return;
   end

   A_out = cell(K + 1, K + 1);
   % Keep variable blocks separate
   A_out(1:K, 1:K) = A(idxMain, idxMain);

   % Condense coupling to saddle interface block (column K+1 and row K+1)
   for i = 1:K
      A_out{i, K + 1} = cell2matrix(A(idxMain(i), idxSupp));
      A_out{K + 1, i} = cell2matrix(A(idxSupp, idxMain(i)));
   end

   % Condense interface-interface block
   A_out{K + 1, K + 1} = cell2matrix(A(idxSupp, idxSupp));
end