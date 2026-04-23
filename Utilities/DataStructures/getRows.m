function out = getRows(A, r)

  if isa(A, 'ArrayOfArrays')
    out = A.getRows(r);
  else
    out = A(r,:);
  end

end

