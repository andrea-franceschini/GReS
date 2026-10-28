function [nEnts, ents, entsPosition] = readEntitySet(input, mesh, entityType, bcName)
% read entity set (can be a scalar tag or a path to a file)

if isfield(input,"bcListFile")
  % input list is in a file
  assert(isscalar(fieldnames(input)), ...
    "Error in BC %s: A path to the entityList should be the unique input of the " + ...
    "BCentities field. ", bcName);
  entPath = getXMLData(input,[],"bcListFile");
  [nEnts, ents] = readListFile(entPath);

  entsPosition = getLocation(ents,mesh,entityType);

  return

else
  entsID = [];
  if isfield(input,"surfaceTags")
    % input is a surface tag with component specification
    surfTags = getXMLData(input,[],"surfaceTags");
    switch entityType
      case "NodeBC"
        entsID = unique(mesh.surfaces(ismember(mesh.surfaceTag,surfTags),:));
      case "SurfBC"
        entsID = find(ismember(mesh.surfaceTag,surfTags));
      otherwise
        error("Error for BC %s: XML field surfaceTags is not valid for BC of type %s", bcName, entityType)
    end
  elseif isfield(input,"bcList")
    % direct entity assignment in the xml file
    entsID = getXMLData(input,[],"bcList");
  elseif isfield(input,"box")
    boxSize = getXMLData(input,[],"box");
    Lx = boxSize(1:2);
    Ly = boxSize(3:4); 
    Lz = boxSize(5:6);
    switch entityType
      case "NodeBC"
        c = mesh.coordinates;
      case {"ElementBC","VolumeForce"}
        c = mesh.cellCentroid;
      otherwise
          error("Error for BC %s: XML field box is not valid for BC of type %s", bcName, entityType)
    end

    entsID = all([ c(:,1) > Lx(1), c(:,1) < Lx(2),...
                   c(:,2) > Ly(1), c(:,2) < Ly(2),...
                   c(:,3) > Lz(1), c(:,3) < Lz(2)],2);
      
    entsID = find(entsID);
  end

  if isempty(entsID)
    error("Error for BC %s: Empty or invalid list of entity in Boundary condition input.",bcName)
  end

  compID = true;

  if isfield(input,"components")
    dir = getXMLData(input,[],"components");
    compID =  ismember(["x","y","z"],dir);
    if ~any(compID)
      error("Error for BC %s: Check syntax of 'components' field.",bcName)
    end
  end

  nEnts = numel(entsID).*compID;
  ents = repmat(entsID,sum(compID),1);
  ents = reshape(ents,[],1);

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
    nEnts = sscanf(line, '%i');
    header = true;
  end
end
if (~header)
  error('Missing header in readSet.');
end
nValMax = sum(nEnts);
ents = zeros(nValMax,1);
id = 1;
while ~feof(fid)
  line = fgetl(fid);
  word = sscanf(line, '%s');
  if (~strcmp(word(1), '%'))
    % If this is not a commented line (not starting with %)
    num = sscanf(line, '%i');
    nNum = length(num);
    ents(id:id+nNum-1) = num;
    id = id + nNum;
  end
end
fclose(fid);

end


function pos = getLocation(ents,mesh,entType)

switch entType
  case "NodeBC"
    pos = mesh.coordinates(ents,:);
  case "SurfBC"
    pos = mesh.surfaceCentroid(ents,:);
  case {"ElementBC","VolumeForce"}
    pos = mesh.cellCentroid(ents,:);
end
end
