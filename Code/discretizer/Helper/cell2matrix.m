function mat = cell2matrix(mat)

mat = processCellMatrix(mat);

% Finally call the built-in cell2mat
mat = cell2mat(mat);


end

