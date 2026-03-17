function inputStruct = readInput(varargin)
%READINPUT Read input as struct, key-value pairs, or XML file.
%
%   readInput(structIn)
%   readInput('k1',v1,'k2',v2)
%   readInput('fileName.xml')
%   readInput(defaultStruct,userInput)
%
% If nargin > 1 and the first argument is a struct, it is used as
% default parameters.
%
% Nested structures are not validated recursively. For nested defaults,
% call readInput again at the lower level.
%
% OUTPUT:
%   inputStruct : (possibly nested) struct

default = [];
usrIn   = varargin;

% Detect default struct
if nargin > 1 && isstruct(varargin{1})
    default = varargin{1};
    usrIn   = varargin(2:end);
end

% Handle empty string input
if isscalar(usrIn) && isstring(usrIn{1}) && strlength(usrIn{1}) == 0
    usrIn = {};
end

% No user input → return default
if isempty(usrIn) || all(cellfun(@isempty, usrIn))
  inputStruct = default;
  return
end

% Interpret user input
if numel(usrIn) > 1
    usr = readKeyValueInput(usrIn);
else
    arg = usrIn{1};

    if isstruct(arg)
        if isfield(arg,'fileName') && ~isempty(fieldnames(arg))
            usr = readXMLfile(arg.fileName);
        else
            usr = arg;
        end
    else
        usr = readXMLfile(arg);
    end
end

% Merge with default (if provided)
if isempty(default)
    inputStruct = usr;
else
    inputStruct = mergeInput(default, usr);
end

end






% HELPER FUNCTIONS

function str = readKeyValueInput(kv)
% Convert key-value cell input into a struct

if mod(numel(kv),2) ~= 0
    error('readKeyValueInput:InvalidInput', ...
        'Key-value input must have even length.');
end

str = struct();

for k = 1:2:numel(kv)
    key = kv{k};
    val = kv{k+1};

    if iscell(val) && mod(numel(val),2) == 0
        val = kv2struct(val);   % assumed external utility
    end

    str.(key) = val;
end

end


function str = readXMLfile(fileName)
% Convert XML file into MATLAB struct

str = readstruct(fileName, AttributeSuffix="");
str = parseXMLStruct(str);

end


function S = parseXMLStruct(S)
% Recursively reinterpret values of a struct read from XML

if ~isstruct(S)
    S = getXMLValue(S);   % assumed external utility
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
% Merge default and user-provided structs (top-level only).
%
% Rules:
% 1. usr overrides default.
% 2. Missing fields are inherited.
% 3. Empty default field → required in usr.
% 4. Enforce type consistency when both are non-empty.

out  = default;
fusr = fieldnames(usr);
fdef = fieldnames(default);

% Apply user fields
for k = 1:numel(fusr)

    fn   = fusr{k};
    uval = usr.(fn);

    if isfield(default, fn)
        dval = default.(fn);

        % Required field check
        if isempty(dval) && isempty(uval) 
            error('mergeInput:RequiredFieldNotProvided', ...
                'Required field "%s" must be provided.', fn);
        end

        % Type consistency
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

% Check required fields missing in usr
for k = 1:numel(fdef)

    fn = fdef{k};

    if isempty(default.(fn)) && ~isfield(usr, fn)
        error('mergeInput:RequiredFieldMissing', ...
            'Required field "%s" must be provided.', fn);
    end
end

end


