function [A] = fixPattern(A)
% FIXPATTERN Ensures block matrix A has non-empty blocks of valid dimensions
% and symmetric nonzero sparsity patterns.

   N = size(A, 1);

   % Determine row and column dimensions for each block row/col index.
   % Also flag blocks that are empty [] or 0x0.
   blkRows = zeros(N, 1);
   blkCols = zeros(N, 1);
   emptyMask = false(N, N);

   for j = 1:N
      for i = 1:N
         M = A{i,j};
         if isempty(M)
            emptyMask(i,j) = true;
         else
            [r, c] = size(M);
            if r > 0 && blkRows(i) == 0, blkRows(i) = r; end
            if c > 0 && blkCols(j) == 0, blkCols(j) = c; end
         end
      end
   end

   % For square block partitionings, if an entire row or column had only empty
   % blocks, recover the missing dimension from the symmetric counterpart.
   for i = 1:N
      if blkRows(i) == 0 && blkCols(i) > 0
         blkRows(i) = blkCols(i);
      elseif blkCols(i) == 0 && blkRows(i) > 0
         blkCols(i) = blkRows(i);
      end
   end

   % Replace empty blocks with correctly dimensioned sparse zeros.
   if any(emptyMask, 'all')
      [eRows, eCols] = find(emptyMask);
      for k = 1:numel(eRows)
         i = eRows(k);
         j = eCols(k);
         r = blkRows(i);
         c = blkCols(j);
         if r == 0 || c == 0
            error('fixPattern:cannotDetermineBlockSize', ...
                  'Cannot determine dimensions for empty block (%d,%d).', i, j);
         end
         A{i,j} = sparse(r, c);
      end
   end

   % Symmetrize the sparsity pattern across blocks.
   % If A{i,j} has a nonzero at (r,c) but A{j,i} does not at (c,r), add an
   % eps-weight entry at (c,r) in A{j,i} so the graph structure is symmetric.
   for j = 1:N
      for i = 1:j
         Aij = A{i,j};
         if i == j
            % Diagonal block: check internal symmetry of sparsity pattern
            if nnz(Aij) > 0
               S = spones(Aij);
               St = S';
               if ~isequal(S, St)
                  % Nonzeros in S not present in S'
                  diff = S - St;
                  [r, c] = find(diff == 1);
                  A{i,i} = Aij + sparse(c, r, eps, blkRows(i), blkRows(i));
               end
            end
         else
            % Off-diagonal block pair (i,j) and (j,i)
            Aji = A{j,i};
            nnz_ij = nnz(Aij);
            nnz_ji = nnz(Aji);
            % Fast exit if both blocks are purely zero
            if nnz_ij == 0 && nnz_ji == 0
               continue;
            end

            Sij = spones(Aij);
            Sji_t = spones(Aji)';

            % Fast equality check: compares sparsity structures directly
            if ~isequal(Sij, Sji_t)
               diff = Sij - Sji_t;
               % Entries present in A{i,j} but missing in A{j,i}'
               [r1, c1] = find(diff == 1);
               if ~isempty(r1)
                  A{j,i} = Aji + sparse(c1, r1, eps, blkCols(j), blkRows(i));
               end
               % Entries present in A{j,i}' but missing in A{i,j}
               [r2, c2] = find(diff == -1);
               if ~isempty(r2)
                  A{i,j} = Aij + sparse(r2, c2, eps, blkRows(i), blkCols(j));
               end
            end
         end
      end
   end
end
