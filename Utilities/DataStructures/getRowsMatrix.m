function out = getRowsMatrix(A, r)
  if isa(A, 'ArrayOfArrays')
    out = A.getRowsMatrix(r);
  else
    out = A(r,:);
  end
end

