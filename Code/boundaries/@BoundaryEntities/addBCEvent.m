function addBCEvent(obj,varargin)

default = struct('time',[],...
  'value',[]);

params = readInput(default,varargin{:});

tVal = params.time;

if ~isnumeric(tVal)
  time = -1;
else
  time = tVal;
end

obj.times(end+1) = time;
obj.bcData = [obj.bcData; struct('time', tVal, 'value', params.value)];

% reorder time in ascending order
[obj.times,s] = sort(obj.times,"ascend");
obj.bcData = obj.bcData(s);

obj.nTimes = sum(obj.times >= 0.0);

% check for repeated bc times
if length(unique(obj.times))~=length(obj.times)
  error("Multiple BC events with same time are" + ...
    "not allowed")
end

end
