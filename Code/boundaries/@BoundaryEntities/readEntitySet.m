function [nEnts, ents, entsPosition] = readEntitySet(obj,type, ents, components, mesh)
% read entity set and return the number of entities for each
% component, the list of entities and their reference location

switch lower(type)
  % input file for list of entities
  case "bclistfile"
    [nEnts, ents] = obj.readListFile(ents);
    entsPosition = getLocation(obj,ents,mesh);
    return
  case {'tags','tag'}
    tags = ents;
    entsID = getEntitiesFromTags(obj.sourceField,mesh,obj.sourceField,tags);
  case "bclist"
    entsID = ents;
  case "box"
    boxSize = ents;
    Lx = boxSize(1:2);
    Ly = boxSize(3:4);
    Lz = boxSize(5:6);
    switch obj.sourceField
      case "node"
        c = mesh.coordinates;
      case "cell"
        c = mesh.cellCentroid;
      otherwise
        error("Error for BC %s: entityListType 'box' is not valid for BC of type %s", obj.name, obj.sourceField)
    end

    entsID = all([ c(:,1) > Lx(1), c(:,1) < Lx(2),...
      c(:,2) > Ly(1), c(:,2) < Ly(2),...
      c(:,3) > Lz(1), c(:,3) < Lz(2)],2);

    entsID = find(entsID);

  otherwise
    error("Unrecognized entityListType '%s' for Boundary condition '%s'.\n" + ...
      "Valid entries for 'entityListType' are: 'bclist','bclistfile','tags','box'.",type,obj.name)
end


if isempty(entsID)
  error("Error for BC %s: Empty or invalid list of entity in Boundary condition input.",obj.name)
end

% expand entity list to components
compID = true;

if ~isempty(components)

  if isnumeric(components)
    dir = ["x","y","z"];
    components = dir(components);
  end

  compID =  ismember(["x","y","z"],components);

end

if ~any(compID)
  error("Error for BC %s: Check syntax of 'components' field.",obj.name)
end

nEnts = numel(entsID).*compID;
entsID = reshape(entsID,[],1);
ents = repmat(entsID,sum(compID),1);
entsPosition = getLocation(obj,ents,mesh);

end

