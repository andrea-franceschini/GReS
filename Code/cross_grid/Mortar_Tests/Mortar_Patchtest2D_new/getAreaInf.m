function [nId,lInf] = getAreaInf(mesh,tag)
% given an edge tag, return:
% nId: node entries in sorted order
% lInf: length influence for each node
nList = unique(mesh.edges(mesh.edgeTag == tag,:));
% get distance from origin as measure for sorting
dist = sqrt(mesh.coordinates(nList,1).^2+mesh.coordinates(nList,2).^2);
[~,id] = sort(dist,'ascend');
diffPos = diff(mesh.coordinates(nList(id),:));
lInf = sqrt(diffPos(:,1).^2+diffPos(:,2).^2+diffPos(:,3).^2);
lInf = 0.5*[lInf(1); lInf(1:end-1)+lInf(2:end); lInf(end)];
nId = nList(id);
end

