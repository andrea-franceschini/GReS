function gpIds = getGaussPointIds(map, cellIds)

% map: a map from cell index to gp index with associated number of gp
% ngp: stores number of gp for each cell
% cellIds: the cells for which expansion is requested

if size(map,2) ~=2
  error("Input map must have the same length as the number of gauss points array")
end

counts = map(cellIds(:),2);
starts = map(cellIds(:),1);

gpIds = repelem(starts, counts,1) ...
      + (0:sum(counts)-1).' ...
      - repelem(cumsum([0; counts(1:end-1)]), counts,1);

end