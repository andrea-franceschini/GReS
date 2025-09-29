function [N, pos] = computeBasisF1D(elem_ID, nInts, topol, nodes, degree)
% Find the value that a linear basis functions take at some points defined in
% a list
% each column store the value of one basis function over each point
% intPts = [-1 -0.99 0.99 1];
% intPts = unique([intPts, linspace(intPts(2), intPts(4), nInts)]);
switch degree
    case 0 % linear interpolation of 2nd order elements for contact detection
        intPts = linspace(-1,1,4);
        N = 0.5 - 0.5*intPts';
        i1 = nodes(topol(elem_ID,1),:);
        i2 = nodes(topol(elem_ID,3),:);
    case 1
        intPts = linspace(-1,1,nInts);
        N = zeros(nInts,2);
        N(:,1) = 0.5 - 0.5*intPts';
        N(:,2) = 0.5 + 0.5*intPts';
        i1 = nodes(topol(elem_ID,1),:);
        i2 = nodes(topol(elem_ID,2),:);
    case 2
        % the interpolation honor also the midpoint ( corresponding to
        % internal nodes)
        intPts = linspace(-1,1,nInts);
        N = computeQuadraticSF(intPts');
        i1 = nodes(topol(elem_ID,1),:);
        i2 = nodes(topol(elem_ID,3),:);
    otherwise
        error('Basis functions of order > 2 are not supported')
end
% get coordinates of interpolation points in the real space
pos = ref2nod(intPts, i1, i2);
end


