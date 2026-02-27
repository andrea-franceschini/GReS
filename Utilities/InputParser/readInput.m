function inputStrcut = readInput(varargin)
% Read input as struct, key-value pair or xml input file.
% if nargin > 1 and the first input is a struct, this is used as default
% set of parameters.
% In case of nested input structure, the default check happens only at the
% higher level. Inner structure checks must be done by calling readInput
% again, with a lower-level set of input parameters
%  EXAMPLE:
%   readInput(structIn)
%   readInput('k1',v1,'k2',v2)
%   readInput('fileName.xml')
%   readInput(defaultStruct,userInput)

% Return a (possibly nested) struct

default = [];

if nargin > 1 && isstruct(varargin{1})
  default = varargin{1};
  usrIn = varargin(2:end);
else
  usrIn = varargin(:);
end

if numel(usrIn) > 1
  usr = readKeyValueInput(usrIn);
else
  if isstruct(usrIn{:})
    usr = usrIn{:};
  else
    usr = readXMLfile(usrIn{:});
  end
end

if isempty(default)
  inputStrcut = usr;
else
  inputStrcut = mergeInput(default,usr);
end

end



% --------------- HELPER FUNCTIONS ----------------

function str = readKeyValueInput(kv)
% convert a key-value input into a structure
if mod(numel(kv),2) ~= 0
  error('Key-value cell must have even length');
end

str = struct();
for k = 1:2:numel(kv)
  key = kv{k};
  val = kv{k+1};

  if iscell(val) && mod(numel(val),2) == 0
    val = kv2struct(val);
  end

  str.(key) = val;
end
end

function str = readXMLfile(fileName)
% convert xml file into a matlab struct

str = readstruct(fileName,AttributeSuffix="");

% recursively interpret XML data
str = parseXMLStruct(str);
end

function S = parseXMLStruct(S)
% Reintepret values of a struct read from xml file

if ~isstruct(S)
  % Interpret and read xml data
  S = getXMLValue(S);
  return;
end

% Struct (possibly array)
for i = 1:numel(S)
  fn = fieldnames(S(i));
  for k = 1:numel(fn)
    field = fn{k};
    S(i).(field) = parseXMLStruct(S(i).(field));
  end
end
end


function out = mergeInput(default, usr)
% Merge default and user-provided structs (shallow, top-level only)
% To merge lower-level struct, rerun the command
%
% Rules:
% 1. Fields in usr override default.
% 2. If a field is missing in one, take the one present.
% 3. If a field in default is [], it is required in usr.
% 4. If both default and usr have values, enforce type consistency.

out = default;
fusr = fieldnames(usr);
fdef = fieldnames(default);

for k = 1:numel(fusr)
  fn = fusr{k};
  uval = usr.(fn);

  if isfield(default, fn) 
    dval = default.(fn);

    % If default is empty, user must provide
    if isempty(dval) && isempty(uval)
      error('mergeInput:RequiredFieldNotProvided', ...
        'Required field "%s" must be provided by the user.', fn);
    end

    % Type check if both are non-empty non-missing
    if ~isempty(dval) && ~isempty(uval) && ~any(ismissing(dval))
      if ~isa(uval, class(dval))
        error('mergeInput:TypeMismatch', ...
          'Field "%s" has type %s in user input but expected %s from default.', ...
          fn, class(uval), class(dval));
      end
    end

  end

  out.(fn) = uval; 
end

% Check required fields from default that were missing in usr
for k = 1:numel(fdef)
  fn = fdef{k};
  if isempty(default.(fn)) && ~isfield(usr, fn)
    error('mergeInput:RequiredFieldMissing', ...
      'Required field "%s" must be provided by the user.', fn);
  end
end
end


