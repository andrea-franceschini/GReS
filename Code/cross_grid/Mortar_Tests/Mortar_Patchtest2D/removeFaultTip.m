function [mesh,nodes] = removeFaultTip(mesh,nodes, N)
% Remove the last N edges near the fault boundary
[~,id] = sort(mesh.coordinates(nodes,1));
idDel = nodes([id(1:N) id(end+1-N:end)]); % list of nodes to be removed
id1 = any(ismember(mesh.edges,idDel),2);
id2 = mesh.edgeTag == 1;
idCellDel = all([id1 id2],2);
mesh.edges = mesh.edges(~idCellDel,:);
mesh.edgeTag = mesh.edgeTag(~idCellDel,:);
nodes = nodes(~ismember(nodes,idDel));
end

