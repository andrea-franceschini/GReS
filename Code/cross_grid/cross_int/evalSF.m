function [vals, pos] = evalSF(nodeID, elem_list, nInts, mastertop, master)
% evaluate shape function in the real space and return position of
% integration points in the real space (extended to x,y)
% already ordered to perform RBF interpolation
intPts = [-1 -0.9999 0.9999 1];
intPts = unique([intPts, linspace(intPts(2), intPts(4), nInts)]);
vals = zeros(length(intPts)*length(elem_list),1);
pos =  zeros(length(intPts)*length(elem_list),2);
k = 0;
for i = elem_list'
    if nodeID == mastertop(i,1)
        % shape function equal to 1 in the right node
        vals(k+1:k+length(intPts)) = 0.5 - 0.5*intPts';
    else
        % shape function equal to 1 in the right node
        vals(k+1:k+length(intPts)) = 0.5 + 0.5*intPts';
    end
    i1 = master(mastertop(i,1),:);
    i2 = master(mastertop(i,2),:);
    % get interpolation points in the real space
    pos(k+1:k+length(intPts),:) = ref2nod(intPts, i1, i2);
    k = k + length(intPts);
end

% need to remove duplicate points on selected node
vals = [vals(1:nInts+2); vals(nInts+4:end)];
pos = [pos(1:nInts+2,:); pos(nInts+4:end,:)];
end

