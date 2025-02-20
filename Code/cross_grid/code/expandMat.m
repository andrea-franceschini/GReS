function mat_new = expandMat(mat, N)
    % Get the size of the original matrix
    [s1, s2] = size(mat);
    
    % Initialize the sparse matrix: row indices, column indices, and values
    rows = [];
    cols = [];
    values = [];
    
    % Loop through the original matrix and populate the sparse matrix
    for s = N-1:-1:0
       % Get the row and column indices for the block
       r1 = N*(1:s1) - s;
       r2 = N*(1:s2) - s;
       [colIdx,rowIdx]  = meshgrid(r2,r1);
       % Add the values from the original matrix to the sparse matrix
       rows = [rows; rowIdx(:)];
       cols = [cols; colIdx(:)];
       values = [values; mat(:)];
    end
    
    % Create the sparse matrix directly from the row, column indices, and values
    mat_new = sparse(rows, cols, values, s1 * N, s2 * N);
end