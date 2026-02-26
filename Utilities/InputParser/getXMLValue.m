function val = getXMLValue(raw)
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

if ismissing(raw)
  val = [];
  return
end

error('getXMLValue:invalidType', ...
  'Could not interpret XML value of type %s', class(raw));
end