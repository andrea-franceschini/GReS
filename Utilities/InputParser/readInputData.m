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


