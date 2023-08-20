function [entities_list,surf_list,loaded_surf] = writeBC(grid,table,bclist)
nbc = length(table);
entities_list = cell(2,nbc);
surf_list = cell(2,nbc);

for i = 1:nbc
    list = [];
    surf = cell2mat(table(i));
    for j = 1:length(surf)
    tmp = grid.topology.surfaces(grid.topology.surfaceTag == surf(j),:);
    tmp = tmp(:);
    list = [list;tmp];
    end
    list = unique(list);
    entities_list(:,i)= {bclist(i);list};
end
for i = 1:nbc
    listSurf = [];
    surf = cell2mat(table(i));
    for j = 1:length(surf)
        tmp = find(grid.topology.surfaceTag == surf(j));
        tmp = tmp(:);
        listSurf = [listSurf;tmp];
    end
    listSurf = unique(listSurf);
    surf_list(:,i)= {bclist(i);listSurf};
end


%provisional code to define surfaces inside loading area (10m radius)
%BC for surface load must be numb. 2;
% surf_top = cell2mat(surf_list(2,2));
% 
% list_coord = grid.topology.coordinates(grid.topology.surfaces(surf_top,:)',1:2);
% list_coord_reshape = reshape(list_coord',2,3,[]);
% list_coord_sum = sum(list_coord_reshape,2)/3;
% list_coord_column = reshape(list_coord_sum,[],1);
% top_centroids = accumarray(repelem(1:length(surf_top),2)',list_coord_column.^2);
% top_centroids = sqrt(top_centroids);
% 
% loaded_surf= surf_top(top_centroids<10.1);








end

