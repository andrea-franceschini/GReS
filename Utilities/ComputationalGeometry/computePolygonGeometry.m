function varargout = computePolygonGeometry(poly, nV, polyNormals)

% mex accelerated computational geometry for general polygons

% normals are not provided, the code assumes that vertices are already
% given in perimetrical order (no crossing edges)

% if normal is provided, vertices can be a generic set of coplanar point
% and proper CCW ordering is done

if nargin == 2
  [areas,centers,normals] = mxPolygonGeometry(poly, nV);
else
  [areas,centers,normals] = mxPolygonGeometry(poly, nV, polyNormals);
end

varargout{1} = areas;
varargout{2} = centers;

if nargout > 2
  varargout{3} = normals;
end

end

