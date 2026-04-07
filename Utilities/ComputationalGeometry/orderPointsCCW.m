function varargout = orderPointsCCW(poly)

perm = mxOrderPointsCCW(poly);
varargout{1} = poly(perm,:);

if nargout > 1
  varargout{2} = perm;
end

end

