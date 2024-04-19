function [vals, pos] = evalSF_3D(nodeID, elem_list, nInts, mastertop, master, elem)
% evaluate shape function in the real space and return position of
% integration points in the real space (extended to x,y)
% already ordered to perform RBF interpolation
intPts = [-0.999 0.999]; % ad-hoc interpolation
intPts = unique([intPts, linspace(intPts(1), intPts(2), nInts)]);
[y, x] = meshgrid(intPts, intPts);
intPts = [x(:), y(:)];

nSupp = unique(mastertop(elem_list,:));
vals = zeros(size(intPts,1)*length(elem_list),1);
pos =  zeros(size(intPts,1)*length(elem_list),3);
k = 0;
% get value of the basis function in the inner of elem of the support
for i = elem_list'
    % compute basis functions at interpolation points (all 4 nodes)
    bf = computeBasisF(elem.quad,intPts);
    % get basis functions of the given node 
    vals(k+1:k+length(intPts),:) = bf(:,mastertop(i,:) == nodeID);
    % get coords of interpolation points in the real space
    pos(k+1:k+length(intPts),:) = bf*master(mastertop(i,:),:);
    k = k + length(intPts);
end
% add values and positions of the nodes of the support
vals = [vals; 1; zeros(length(nSupp)-1,1)];
posSupp = [master(nSupp(nSupp==nodeID),:);
    master(nSupp(nSupp~=nodeID),:)];
pos = [pos; posSupp];
end

