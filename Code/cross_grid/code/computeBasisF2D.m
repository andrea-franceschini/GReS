function [vals, pos] = computeBasisF2D(elemID, nInts, topol, coord, elem)
% evaluate shape function in the real space and return position of
% integration points in the real space (extended to x,y)
% already ordered to perform RBF interpolation
intPts = [-1 1]; % ad-hoc interpolation
intPts = linspace(intPts(1), intPts(2), nInts);
[y, x] = meshgrid(intPts, intPts);
intPts = [x(:), y(:)];

% get value of the basis function in the inner of elem of the support
% compute basis functions at interpolation points (all 4 nodes)
bf = computeBasisF(elem.quad,intPts);
% get basis functions of the given node
vals = bf;
% get coords of interpolation points in the real space
pos = bf*coord(topol(elemID,:),:);

% add interpolation points of the nodes (otherwise I'll generate useless
% nodes on the edges
% vals = [vals; eye(4)];
% pos = [pos; coord(topol(elemID,:)',:)];
end

