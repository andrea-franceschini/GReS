function vals = readDataSetNew(bc,varargin)

nEnts = bc.nEntities;
totEnts = bc.totEnts;

if isempty(varargin)
  assert(isscalar(bc.bcData),"The id of the boundary data needs to be specified")
end


% get array of BC values depending on the data

if isnumeric(data) && isscalar(data)
  % scalar value 
  vals = repmat(data,nVals,1);
elseif isValidFunction(data)
  % function of spacial coordinates (constant 

else 
  % file with BC values
  vals = readValuesList(data);
end
end


function tf = isValidFunction(expr)
    tf = true;
    try
        str2func(['@(x,y,z)', expr]);
    catch
        tf = false;
    end
end

function vals = readValuesList(fileName)
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
if length(vals) ~= nVals
  error('Number of values in %s not matching number of constrained entities',fileName)
end
fclose(fid);
end
