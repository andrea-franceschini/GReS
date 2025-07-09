function doflist = getDoF(nodeList)
% get dof numbering corresponding to a list of nodes for 2D problem

doflist = 2*repelem(nodeList,2) + repmat([-1 0],1,length(nodeList));
end

