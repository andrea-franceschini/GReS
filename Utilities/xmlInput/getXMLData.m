function data = getXMLData(str, default, varargin)
% getXMLData - recursively extract XML field from structure
%
% Usage:
%   data = getXMLData(str, default, 'root', 'child', 'subchild', ...)
%
% If any field is missing or empty/"Default", returns the given default value.
% If default is empty, the XML data is required.

if isempty(varargin)
   data = str;
   return;
end

% Current field name
field = varargin{1};

% Check if field exists
if isstruct(str) && isfield(str, field)
   value = str.(field);
   % Recursive call with remaining fields
   data = getXMLData(value, default, varargin{2:end});
else
   data = default;
   if isempty(data)
      error('Input parameter %s is required', field)
   end
   return;
end

% If final data is empty or "Default", replace with default value
if ischar(data) && (strcmp(data, "") || strcmpi(data, "Default"))
   data = default;
elseif isempty(data)
   data = default;
end

if isempty(data)
  error('Input parameter %s is required', field)
end
end