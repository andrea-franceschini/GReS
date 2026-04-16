function inputStruct = readInput(varargin)
%READINPUT Parse input arguments into a struct, with optional defaults.
%
%   inputStruct = readInput(structIn)
%   inputStruct = readInput('k1',v1,'k2',v2,...)
%   inputStruct = readInput('fileName.xml')
%   inputStruct = readInput(defaultStruct, userInput)
%
%   Accepts four input forms:
%     - A single struct, returned as-is (or used to look up an XML file
%       if it contains a 'fileName' field).
%     - Name-value pairs, converted to a struct.
%     - A path to an XML file, parsed into a struct.
%     - A default struct followed by any of the above: user values are
%       merged into the default (user overrides default; see mergeInput).
%     Merge rules:
%     - User values override matching default fields.
%     - Default fields absent in usr are inherited unchanged.
%     - A default field that is empty marks that field as required:
%       usr must supply a non-empty value, otherwise an error is raised.
%     - When both values are non-empty, their classes must match.
%
%   Nested structures are not validated recursively. For nested defaults,
%   call readInput again at the lower level.
%
%   INPUT:
%     varargin : default struct (optional, must be first) followed by one
%                of: struct, name-value pairs, or XML file path.
%
%   OUTPUT:
%     inputStruct : scalar struct with the parsed (and merged) fields.
%
%   See also: mergeInput, readKeyValueInput, readXMLfile.
 
default = [];
usrIn   = varargin;
 
% If first arg is a struct and more args follow, treat it as defaults
if nargin > 1 && isstruct(varargin{1})
    default = varargin{1};
    usrIn   = varargin(2:end);
end
 
% Treat empty string as no input
if isscalar(usrIn) && isstring(usrIn{1}) && strlength(usrIn{1}) == 0
    usrIn = {};
end
 
% No user input, return default unchanged
if isempty(usrIn) || all(cellfun(@isempty, usrIn))
  inputStruct = default;
  return
end

% Parse user input
if numel(usrIn) > 1
  % Multiple args, name-value pairs
  usr = readKeyValueInput(usrIn);
else
  arg = usrIn{1};

  if isstruct(arg)
    if isfield(arg,'fileName') && ~isempty(fieldnames(arg))
      usr = readXMLfile(arg.fileName);   % struct with fileName → XML
    else
      usr = arg;                         % plain struct → pass through
    end
  elseif isstring(arg) || ischar(arg)
    usr = readXMLfile(arg);                % string → XML file path
  else
    % argument is type missing
    usr = [];
  end
end

% Merge with defaults if provided
if isempty(default)
  inputStruct = usr;
else
    inputStruct = mergeInput(default, usr);
end
 
end
 



function str = readKeyValueInput(kv)
%READKEYVALUEINPUT Convert a name-value cell array to a struct.
%
%   str = readKeyValueInput(kv)
%
%   kv must have even length: {key1, val1, key2, val2, ...}.
%   If a value is itself an even-length cell array, it is recursively
%   converted via kv2struct (external utility).
%
%   INPUT:
%     kv  : cell array of alternating string keys and values.
%
%   OUTPUT:
%     str : scalar struct with one field per key.
 
if mod(numel(kv),2) ~= 0
    error('readKeyValueInput:InvalidInput', ...
        'Key-value input must have even length.');
end
 
str = struct();
 
for k = 1:2:numel(kv)
    key = kv{k};
    val = kv{k+1};
 
    % if iscell(val) && mod(numel(val),2) == 0
    %     val = kv2struct(val);   % recurse into nested name-value cell
    % end
 
    str.(key) = val;
end
 
end
 
 
function str = readXMLfile(fileName)
%READXMLFILE Read an XML file and return its contents as a struct.
%
%   str = readXMLfile(fileName)
%
%   Uses readstruct with no attribute suffix, then post-processes values
%   via parseXMLStruct to convert XML-typed scalars to native MATLAB types.
%
%   INPUT:
%     fileName : path to XML file (char or string).
%
%   OUTPUT:
%     str : struct mirroring the XML element hierarchy.
 
str = readstruct(fileName, AttributeSuffix="");
str = parseXMLStruct(str);
 
end
 
 
function S = parseXMLStruct(S)
%PARSEXMLSTRUCT Recursively convert XML-typed values to MATLAB types.
%
%   S = parseXMLStruct(S)
%
%   Walks every field of S depth-first. Non-struct leaves are converted
%   via getXMLValue (external utility). Struct arrays are handled
%   element-by-element.
%
%   INPUT / OUTPUT:
%     S : struct (possibly nested/array) or scalar XML value.
 
if ~isstruct(S)
    S = getXMLValue(S);   % convert leaf value (external utility)
    return
end
 
for i = 1:numel(S)
    fn = fieldnames(S(i));
    for k = 1:numel(fn)
        field = fn{k};
        S(i).(field) = parseXMLStruct(S(i).(field));
    end
end
 
end
 
 
function out = mergeInput(default, usr)
%MERGEINPUT Merge user struct into a default struct (top-level fields only).
%
%   out = mergeInput(default, usr)
%
%   Merge rules:
%     - User values override matching default fields.
%     - Default fields absent in usr are inherited unchanged.
%     - A default field that is empty marks that field as required:
%       usr must supply a non-empty value, otherwise an error is raised.
%     - When both values are non-empty, their classes must match.
%
%   INPUT:
%     default : struct with baseline field values ([] marks required).
%     usr     : struct with user-supplied overrides.
%
%   OUTPUT:
%     out     : merged struct.
 
out  = default;

if any([isempty(usr),any(ismissing(usr))])
    return
end
fusr = fieldnames(usr);
fdef = fieldnames(default);
 
% Apply user-supplied fields
for k = 1:numel(fusr)
 
    fn   = fusr{k};
    uval = usr.(fn);
 
    if isfield(default, fn)
        dval = default.(fn);
 
        % Required field: default is empty, user must provide a value
        if isempty(dval) && isempty(uval)
            error('mergeInput:RequiredFieldNotProvided', ...
                'Required field "%s" must be provided.', fn);
        end
 
        % Type check: classes must agree when both values are non-empty
        if ~isempty(dval) && ~isempty(uval) && ~any(ismissing(dval))
            if ~isa(uval, class(dval))
                error('mergeInput:TypeMismatch', ...
                    'Field "%s" has type %s but expected %s.', ...
                    fn, class(uval), class(dval));
            end
        end
    end
 
    out.(fn) = uval;
end
 
% Check that all required default fields (empty) were supplied by usr
for k = 1:numel(fdef)
 
    fn = fdef{k};
 
    if isempty(default.(fn)) && ~isfield(usr, fn)
        error('mergeInput:RequiredFieldMissing', ...
            'Required field "%s" must be provided.', fn);
    end
end
 
end