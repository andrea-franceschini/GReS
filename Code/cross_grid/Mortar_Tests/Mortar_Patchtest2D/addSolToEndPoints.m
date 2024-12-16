function u = addSolToEndPoints(mesh,nodes,dofM,dofS,dofIm,dofIs,u)
[~,id] = sort(mesh.coordinates(nodes,1),'ascend');
n1 = nodes(id(1));
n2 = nodes(id(2));
n3 = nodes(id(end-1));
n4 = nodes(id(end));
tmp = [n1 n2; n4 n3];
% sum rows and column of end points to adjacent contribution
for i = 1:2
    nEnd = tmp(i,1);
    nAdj = tmp(i,2);
    for j = -1:0
        % get corresponding dof
        dofEnd = getGlobalDofs(2*nEnd+j,dofM,dofS,dofIm,dofIs,'lagrange');
        dofAdj = getGlobalDofs(2*nAdj+j,dofM,dofS,dofIm,dofIs,'lagrange');
        u(dofEnd) = u(dofAdj);
    end
end
end

