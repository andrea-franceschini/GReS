function mat = processCellMatrix(mat)

% processCellMatrix: makes a sparse block cell array compatible with cell2mat.
%
% - Removes rows/columns that contain only empty cells.
% - Replaces empty cells in remaining rows/columns with sparse matrices
%   of the correct size for their block row and block column.

[nR, nC] = size(mat);

% detect empty row and columns
hasRow = false(nR,1);
hasCol = false(1,nC);

for i = 1:nR
  for j = 1:nC
    if ~isempty(mat{i,j})
      hasRow(i) = true;
      hasCol(j) = true;
    end
  end
end

% remove empty rows/cols
mat = mat(hasRow,hasCol);
[nR,nC] = size(mat);

% determine block size
rowSize = zeros(nR,1);
for i = 1:nR
  for j = 1:nC
    if ~isempty(mat{i,j})
      [rowSize(i),~] = size(mat{i,j});
      break
    end
  end
end

% block column sizes (number of cols in each block column)
colSize = zeros(1,nC);
for j = 1:nC
  for i = 1:nR
    if ~isempty(mat{i,j})
      [~,colSize(j)] = size(mat{i,j});
      break
    end
  end
end

% replace empty blocks with empty sparse matrix
for i = 1:nR
  for j = 1:nC
    if isempty(mat{i,j})
      mat{i,j} = sparse(rowSize(i), colSize(j));
    end
  end
end

end


