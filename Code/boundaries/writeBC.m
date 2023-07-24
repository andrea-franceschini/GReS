function [entities_list] = writeBC(grid,table,bclist)
nbc = length(table);
entities_list = cell(2,nbc);

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

