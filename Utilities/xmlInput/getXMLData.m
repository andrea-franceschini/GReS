function data = getXMLData(fieldStruct, defaultValue, dataName)
% getXMLData - read a parameter from an input file
%
% Input:
%
% str is the structure containing the paramters of a certain xml field


% Check if field exists
if ~isfield(fieldStruct, dataName)
   data = defaultValue;
elseif strcmp(getXMLValue([fieldStruct.(dataName)]),"Default")
   data = defaultValue;
else 
  data = getXMLValue(fieldStruct.(dataName));
end

if isempty(data) 
  error('Input parameter %s is required', dataName)
end
end



function val = getXMLValue(raw)
% getXMLValue - Convert raw XML value into appropriate MATLAB type
%
% Converts strings like "1 2 3" -> [1 2 3], 
% "text" -> "text", 
% "a b c" -> ["a" "b" "c"]

  % If already numeric, return directly
  if isnumeric(raw)
    val = raw;
    return;
  end

  % Normalize to string
  if iscell(raw)
    raw = string(raw);
  elseif ischar(raw)
    raw = string(raw);
  end

  % For string scalars
  if isstring(raw) && isscalar(raw)
    str = strtrim(raw);

    % Try numeric
    num = str2num(str); %#ok<ST2NM>
    if ~isempty(num)
      val = num;
      return;
    end

    % Try string array (split on spaces, commas, semicolons)
    if contains(str, [",", ";", " "])
      parts = split(str, {",", ";", " "});
      parts = parts(parts ~= ""); % remove empties
      val = parts;
      return;
    end

    % Otherwise single string
    val = str;
    return;
  end

  % Already a string array
  if isstring(raw)
    val = raw;
    return;
  end

  error('getXMLValue:invalidType', ...
        'Could not interpret XML value of type %s', class(raw));
end