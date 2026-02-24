function data = readInputData(fieldStruct, defaultValue, dataName)
% getXMLData  Retrieve and convert a parameter from an XML-derived structure
%
%   data = getXMLData(fieldStruct, defaultValue, dataName)
%
%   Retrieves the value associated with the XML field `dataName` within the
%   structure `fieldStruct`. If the field is missing or explicitly set to
%   "Default", the value `defaultValue` is returned. Otherwise, the raw XML
%   value is converted to an appropriate MATLAB type using getXMLValue.
%
%   Inputs:
%     fieldStruct   - Matlab Struct containing parsed XML fields.
%     defaultValue  - Fallback value used when the field is absent or marked
%                     as default.
%     dataName      - Name of the field to extract from the struct.
%
%   Output:
%     data          - Parsed parameter value. Raises an error if the resulting
%                     value is empty.
%
%   See also: getXMLValue


% Check if field exists
if ~isfield(fieldStruct, dataName)
   data = defaultValue;
elseif strcmp(getValue([fieldStruct.(dataName)]),"Default")
   data = defaultValue;
else 
  data = getValue(fieldStruct.(dataName));
end

if isempty(data) 
  error('Input parameter %s is required', dataName)
end
end



function val = getValue(raw)
% getXMLValue  Convert a raw XML field value into a MATLAB type
%
%   val = getXMLValue(raw)
%
%   Interprets XML text or numeric content and converts it to a suitable
%   MATLAB representation:
%     - numeric strings → numeric scalars/vectors (e.g. "1 2 3" → [1 2 3])
%     - plain text      → string scalar
%     - space/comma/semicolon-separated tokens → string array
%     - numeric input   → returned unchanged
%
%   Inputs:
%     raw  - Raw XML field content (numeric, char, string, cellstr).
%
%   Output:
%     val  - Interpreted MATLAB value.
%
%   Errors:
%     Throws getXMLValue:invalidType if the input cannot be interpreted.

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
    num = str2num(str);
    if ~isempty(num)
      val = num;
      return;
    end

    % Try string array (split on spaces, commas, semicolons)
    if contains(str, [",", ";", " "])
      parts = regexp(str, '[,; ]+', 'split');
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