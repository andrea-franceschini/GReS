function data = getXMLData(fieldStruct, defaultValue, dataName)
% getXMLData - read a parameter from an input file
%
% Input:
%
% str is the structure containing the paramters of a certain xml field


% Check if field exists
if ~isfield(fieldStruct, dataName)
   data = defaultValue;
elseif strcmp(fieldStruct.(dataName),"Default")
   data = defaultValue;
else 
  data = fieldStruct.(dataName);
end

if isempty(data)
  error('Input parameter %s is required', dataName)
end
end