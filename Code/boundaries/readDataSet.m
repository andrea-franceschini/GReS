function vals = readDataSet(bc,varargin)

%nEnts = bc.nEntities;
%totEnts = bc.totEnts;

if isempty(varargin)
  assert(isscalar(bc.bcData),"Error in BC %s The id of the boundary data needs to be specified",bc.name)
  bcVal = bc.bcData.value;
else
  bcVal = bc.bcData(varargin{1}).value;
end


if isnumeric(bcVal)
  % scalar value
  if isscalar(bcVal)
    vals = repmat(bcVal,bc.totEnts,1);
  else
    vals = bcVal;
  end
elseif isValidFunction(bcVal)
  % function of spatial coordinates
  f = str2func(['@(x,y,z)', char(bcVal)]);
  vals = f(bc.entityPos(:,1),bc.entityPos(:,2),bc.entityPos(:,3));
else 
  % file with BC values
  vals = readValuesList(bcVal,bc.totEnts);
end

assert(bc.totEnts == numel(vals),...
  "Number of BC values not matching number of BC entities")
end


function tf = isValidFunction(expr)
    tf = true;
    try
        f = str2func(['@(x,y,z)', char(expr)]);
        f(1,1,1);
    catch
        tf = false;
    end
end

function vals = readValuesList(fileName,nVals)
% path to input file of numeric values
if (~exist(fileName, 'file'))
  error('File %s does not seem to exist. Please, check the provided file.', fileName);
end
fid = fopen(fileName, 'r');
vals = zeros(nVals,1);
id = 1;
while ~feof(fid)
  line = fgetl(fid);
  word = sscanf(line, '%s');
  if (~strcmp(word(1), '%'))
    % If this is not a commented line (not starting with %)
    num = sscanf(line, '%e');
    nNum = length(num);
    vals(id:id+nNum-1) = num;
    id = id + nNum;
  end
end
vals = vals(1:id-1);
if length(vals) ~= nVals
  error(['Number of values in %s not matching number of constrained' ...
    'degrees of freedom'],fileName)
end
fclose(fid);
end
