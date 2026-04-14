function varargout = orderPointsCCW(poly,normal)

if nargin == 1
  perm = mxOrderPointsCCW(poly);
else
  perm = mxOrderPointsCCW(poly,normal);
end
varargout{1} = poly(perm,:);

if nargout > 1
  varargout{2} = perm;
end

end

