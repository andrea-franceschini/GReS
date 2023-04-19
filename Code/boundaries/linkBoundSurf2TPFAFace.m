function linkBoundSurf2TPFAFace(model,bound,grid)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
keys = bound.db.keys;
flRenum = false(length(keys),1);
for i = 1 : length(keys)
  if strcmp(bound.getCond(keys{i}), 'SurfBC') && model.isFVTPFABased(bound.getPhysics(keys{i}))
    flRenum(i) = true;
    break
  end
end
%
if any(flRenum)
  % Create the face topology as a matrix starting from the vector
  % Select only the faces laying on the domain boundary and sort by row
  %
  faceTop = reshape(grid.faces.nodes2Faces,max(diff(grid.faces.mapN2F)),grid.faces.nFaces);
  ptrBFace = sum(grid.faces.faceNeighbors ~= 0,2) == 1;
  idBFace = find(ptrBFace);
  faceBound = faceTop(:,ptrBFace);
  faceBound = faceBound';
  faceBound = sort(faceBound,2);
end
%
for i = find(flRenum)
  surfID = bound.getDofs(keys{i});
  bFaceTop = grid.topology.surfaces(surfID,:);
  [~,~,ib] = intersect(sort(bFaceTop,2),faceBound,'stable','rows');
  bound.setDofs(keys{i},idBFace(ib));
end

end