% function inputStruct = readXMLfile(fileName)
% % override the MATLAB built-in readstruct to include external xml in it
% 
% inputStruct = readstruct(fileName,AttributeSuffix="");
% 
% % recursively look for <Included> fields in the file and read
% % the included xml file
% inputStruct = processIncludedFields(inputStruct);
% 
% end
% 
% function str = mergeStruct(str,includedStruct)
% flds = fieldnames(includedStruct);
% 
% for i = 1:numel(flds)
%   if isfield(str,flds{i})
%     str.(flds{i}) = includedStruct.(flds{i});
%   else
%     str.(flds{i}) = [str.(flds{i}); includedStruct.(flds{i})];
%   end
% end
% end
% 
% function s = processIncludedFields(s)
% 
% if numel(s) > 1
%   % avoid processing structure array. An Included field must replace a
%   % stand-alone field
%   return
% end
% fn = fieldnames(s);
% 
% for k = 1:numel(fn)
%   field = fn{k};
%   value = s.(field);
% 
%   if isstruct(value)
% 
%     % If this field is 'Included'
%     if strcmp(field, 'Included')
%       nIncluded = numel(value);
% 
%       for iInc = 1:nIncluded
%         incStruct = value(iInc);
% 
%         % Process each file inside Included
%         if isfield(incStruct, 'File')
%           nFile = numel(incStruct.File);
%           for iFile = 1:nFile
%             fileName = incStruct.File(iFile).name;
%             includedData = readXMLfile(fileName);
%             s = mergeStruct(s, includedData);  % update s
%           end
%         end
%       end
%     end
% 
%     % Recurse into nested structs
%     s.(field) = processIncludedFields(value);
%   end
% end
% end


