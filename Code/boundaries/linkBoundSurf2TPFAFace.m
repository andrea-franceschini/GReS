function linkBoundSurf2TPFAFace(model,bound,grid)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
keys = bound.db.keys;
flRenum = false(length(keys),1);
for i = 1 : length(keys)
  if strcmp(bound.getCond(keys{i}), 'SurfBC') && model.isFVTPFABased(bound.getPhysics(keys{i}))
    flRenum(i) = true;
  end
end
%
if any(flRenum)
  % Create the face topology as a matrix starting from the vector
  % Select only the faces lying on the domain boundary and sort by row
  % The RESHAPE functions works only if all faces have the same number of
  % nodes
  %
  faceTop = reshape(grid.faces.nodes2Faces,max(diff(grid.faces.mapN2F)),grid.faces.nFaces);
%   ptrBFace = sum(grid.faces.faceNeighbors ~= 0,2) == 1;
%   idBFace = find(ptrBFace);
%   faceBound = faceTop(:,ptrBFace);
  idBFace = find(sum(grid.faces.faceNeighbors ~= 0,2) == 1);
  faceBound = faceTop(:,idBFace);
  faceBound = sort(faceBound',2);
  %
  for i = find(flRenum)'
    bFaceTop = grid.topology.surfaces(bound.getDofs(keys{i}),:);
    [~,~,ib] = intersect(sort(bFaceTop,2),faceBound,'stable','rows');
    newID = idBFace(ib);
    assert(all(sum(grid.faces.faceNeighbors(newID,:) ~= 0,2) == 1),'Corrupted face renumbering for %s surface condition',bound.getName(keys{i}));
    bound.setDofs(keys{i},idBFace(ib));
  end
end
end