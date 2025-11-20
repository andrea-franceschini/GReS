function linkBoundSurf2TPFAFace(solver)
% TO DO: call this from within the SinglePhaseFlow solver!

bcs = solver.bcs;

keys = bcs.db.keys;
flRenum = false(length(keys),1);
for i = 1 : length(keys)
  if strcmp(bcs.getCond(keys{i}), 'SurfBC') && strcmp(bcs.getVariable(keys{i}),solver.getField())
    flRenum(i) = true;
  end
end

% input solver requires BC renumbering
if any(flRenum)
  % Create the face topology as a matrix starting from the vector
  % Select only the faces lying on the domain boundary and sort by row
  % The RESHAPE functions works only if all faces have the same number of
  % nodes
  %
  faceTop = reshape(solver.faces.nodes2Faces,max(diff(solver.faces.mapN2F)),solver.faces.nFaces);
  idBFace = find(sum(solver.faces.faceNeighbors ~= 0,2) == 1);
  faceBound = faceTop(:,idBFace);
  faceBound = sort(faceBound',2);
  %
  for i = find(flRenum)'
    bFaceTop = solver.mesh.surfaces(bcs.getEntities(keys{i}),:); 
    [~,~,ib] = intersect(sort(bFaceTop,2),faceBound,'stable','rows');
    newID = idBFace(ib);
    assert(all(sum(solver.faces.faceNeighbors(newID,:) ~= 0,2) == 1),'Corrupted face renumbering for %s surface condition',bcs.getName(keys{i}));
    bcs.setDofs(keys{i},idBFace(ib));
  end
end
end