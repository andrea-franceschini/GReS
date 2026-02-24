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
  usrIn = varargin{2:end};
else
  usrIn = varargin{:};
end

if numel(usrIn) > 1
  usr = readKeyValueInput(usrIn);
else
  if isstruct(usrIn{1})
    usr = usrIn{1};
  else
    usr = readXMLfile(usrIn{1});
  end
end

if isempty(default)
  inputStrcut = usr;
else
  inputStrcut = mergeInput(default,usr);
end

end

function str = readKeyValueInput(kv)
if mod(numel(kv),2) ~= 0
  error('Key-value cell must have even length');
end

str = struct();
for k = 1:2:numel(kv)
  key = kv{k};
  val = kv{k+1};

  % Recursively convert even-length cell values into struct
  if iscell(val) && mod(numel(val),2) == 0
    val = kv2struct(val);
  end

  str.(key) = val;
end
end

function str = readXMLfile(fileName)
str = readstruct(fileName);
end


function mergeInput(default,usr)
end


