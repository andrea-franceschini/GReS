function [N, pos] = computeBasisF1D(elem_ID, nInts, topol, nodes, degree)
% Find the value that a linear basis functions take at some points defined in
% a list
% each column store the value of one basis function over each point
intPts = [-1 -0.99 0.99 1];
intPts = unique([intPts, linspace(intPts(2), intPts(4), nInts)]);
switch degree
    case 0 % linear interpolation of 2nd order elements for contact detection
        N = zeros(nInts+2,2);
        N(:,1) = 0.5 - 0.5*intPts';
        N(:,2) = 0.5 + 0.5*intPts';
        i1 = nodes(topol(elem_ID,1),:);
        i2 = nodes(topol(elem_ID,3),:);
    case 1
        N = zeros(nInts+2,2);
        N(:,1) = 0.5 - 0.5*intPts';
        N(:,2) = 0.5 + 0.5*intPts';
        i1 = nodes(topol(elem_ID,1),:);
        i2 = nodes(topol(elem_ID,2),:);
    case 2
        N = computeQuadraticSF(intPts');
        i1 = nodes(topol(elem_ID,1),:);
        i2 = nodes(topol(elem_ID,3),:);
    otherwise
        error('Basis functions of order > 2 are not supported')
end
% get coordinates of interpolation points in the real space
pos = ref2nod(intPts, i1, i2);
end


