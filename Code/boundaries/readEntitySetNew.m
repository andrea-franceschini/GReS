function [nEnts, ents, entsPosition] = readEntitySetNew(input,mesh, entityType)
% read entity set (can be a scalar tag or a path to a file)

if isfield(input,"entityList")
  % input list is in a file
  assert(isscalar(fieldnames(input)), ...
    "A path to the entityList should be the unique input of the " + ...
    "BCentities field. ");
  entPath = getXMLData(input,[],"entityList");
  [nEnts, ents] = readListFile(entPath);

  entsPosition = getLocation(ents,mesh,entityType);

  return
  
else
  if isfield(input,"surfaceTags")
  % input is a surface tag with component specification
  surfTags = getXMLData(input,[],"surfaceTags");
  switch entityType
    case "Nodes"
      entsID = unique(mesh.surfaces(ismember(mesh.surfaceTag,surfTags),:));
    case "Surfaces"
      entsID = find(ismember(mesh.surfaceTag,surfTags));
    otherwise
      error("XML field surfaceTags is not valid for BC of type %s", entityType)
  end
  elseif isfield(input,"entities")
    % direct entity assignment in the xml file
    entsID = getXMLData(input,[],"entities");
  end

  if isempty(entsID)
    error("Invalid list of entity in Boundary condition input")
  end

  compID = true;

  if isfield(input,"components")
    dir = getXMLData(input,[],"components");
    compID =  ismember(["x","y","z"],dir);
  end

  nEnts = numel(entsID).*compID;
  ents = repmat(entsID,sum(compID),1);

  entsPosition = getLocation(ents,mesh,entityType);

end
end

function [nEnts, ents] = readListFile(fileName)

if (~exist(fileName, 'file'))
  error('File %s does not seem to exist. Please, check the provided file.', fileName);
end
header = false;
fid = fopen(fileName, 'r');
while (~feof(fid) && ~header)
  line = fgetl(fid);
  word = sscanf(line, '%s');
  if (~strcmp(word(1), '%'))
    % If this is not a commented line (not starting with %)
    nVals = sscanf(line, '%i');
    header = true;
  end
end
if (~header)
  error('Missing header in readSet.');
end
nValMax = sum(nVals);
vals = zeros(nValMax,1);
id = 1;
while ~feof(fid)
  line = fgetl(fid);
  word = sscanf(line, '%s');
  if (~strcmp(word(1), '%'))
    % If this is not a commented line (not starting with %)
    num = sscanf(line, '%i');
    nNum = length(num);
    vals(id:id+nNum-1) = num;
    id = id + nNum;
  end
end
fclose(fid);

end


function pos = getLocation(ents,mesh,entType)

switch entType
  case "Nodes"
    pos = mesh.coordinates(ents,:);
  case "Surfaces"
    pos = mesh.surfaceCentroid(ents,:);
  case {"Elements","VolumeForce"}
    pos = mesh.cellCentroid(ents,:);
end
end
