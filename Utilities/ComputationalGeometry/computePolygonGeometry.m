function varargout = computePolygonGeometry(poly, nV)

% mex accelerated computational geometry for general polygons

[areas,centers,normals] = mxPolygonGeometry(poly, nV);

varargout{1} = areas;
varargout{2} = centers;

if nargout > 2
  varargout{3} = normals;
end

end

