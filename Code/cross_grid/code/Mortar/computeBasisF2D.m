function [bf,pos,varargout] = computeBasisF2D(elemID, nInts, topol, coord, elem)
% evaluate shape function in the real space and return position of
% integration points in the real space (extended to x,y)
% already ordered to perform RBF interpolation
intPts = [-1 1]; % ad-hoc interpolation
intPts = linspace(intPts(1), intPts(2), nInts);
[y, x] = meshgrid(intPts, intPts);
intPts = [x(:), y(:)];

bf = computeBasisF(elem.quad,intPts);
% get coords of interpolation points in the real space
pos = bf*coord(topol(elemID,:),:);

% get basis functions of lower order element for support detection (only
% for quadrilateral elements)
if nargout==3
    assert(obj.degree==2,'Incorrect number of outputs for basis functions of degree 1')
    varargout = computeBasisF(elem.quadL,intPts);
end
end

