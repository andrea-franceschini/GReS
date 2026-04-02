function vals = getValues(obj, t)
% Computes boundary condition values at time t by interpolating between
% two bracketing data sets. Uses binary search to locate the interval:
%   availStep(1) ------- t --------- availStep(2)
% Avoids redundant reads by caching the last two loaded data sets.

if obj.nTimes == 0
  vals = getValuesScaled(obj, t);
elseif obj.nTimes == 1
  vals = getValuesConstant(obj);
else
  vals = getValuesInterpolated(obj, t);
end

end

% -------------------------------------------------------------------------
function vals = getValuesScaled(obj, t)
% Time-dependent scaling via function handle: val = f(t) * bcVal

if obj.readSetFlag
  obj.availVals(:,1) = readDataSet(obj);
  obj.readSetFlag    = false;
end

tScale = str2func(['@(t)', char(obj.bcData.time)]);
vals   = tScale(t) * obj.availVals(:,1);

end

% -------------------------------------------------------------------------
function vals = getValuesConstant(obj)
% Constant BC over time — reads the data set only once

if obj.readSetFlag
  obj.availVals(:,1) = readDataSet(obj);
  obj.readSetFlag    = false;
end

vals = obj.availVals(:,1);

end


function vals = getValuesInterpolated(obj, t)
% BC with multiple time steps — interpolates between bracketing sets

[i1, i2] = bin_search(obj, t);

if i1 == i2
  vals = getValuesAtEdge(obj, i1);
else
  vals = interpolateBetweenSets(obj, t, i1, i2);
end

end


function vals = getValuesAtEdge(obj, idx)
% Handles the edge case where binary search returns the same index twice

if idx == 1
  obj.availSteps(1) = idx;
  obj.availVals(:,1) = readDataSet(obj, 1);
  vals = obj.availVals(:,1);
elseif idx == obj.nTimes
  obj.availVals(:,2) = readDataSet(obj, obj.nTimes);
  vals = obj.availVals(:,2);
end

end


function vals = interpolateBetweenSets(obj, t, i1, i2)
% Loads the two bracketing data sets (reusing cached slots when possible)
% and linearly interpolates to time t

[p1, p2] = resolveCacheSlots(obj, i1, i2);

fac  = (t - obj.times(i1)) / (obj.times(i2) - obj.times(i1));
vals = obj.availVals(:,p1) + fac * (obj.availVals(:,p2) - obj.availVals(:,p1));

end


function [p1, p2] = resolveCacheSlots(obj, i1, i2)
% Determines which cache slot (1 or 2) holds each bracketing set,
% reading from disk only when a set is not already cached

p1 = find(obj.availSteps == i1);
p2 = find(obj.availSteps == i2);

if isempty(p1) && isempty(p2)           % neither set is cached
  p1 = 1;
  p2 = 2;
  obj.availVals(:,p1) = readDataSet(obj, i1);
  obj.availVals(:,p2) = readDataSet(obj, i2);
  obj.availSteps(p1)  = i1;
  obj.availSteps(p2)  = i2;

elseif ~isempty(p1) && isempty(p2)      % i1 cached, load i2 into the other slot
  p2 = 3 - p1;
  obj.availVals(:,p2) = readDataSet(obj, i2);
  obj.availSteps(p2)  = i2;

elseif ~isempty(p2) && isempty(p1)      % i2 cached, load i1 into the other slot
  p1 = 3 - p2;
  obj.availVals(:,p1) = readDataSet(obj, i1);
  obj.availSteps(p1)  = i1;
end
% if both are already cached, p1/p2 are correct and no reads are needed

end



function [i1, i2] = bin_search(obj, t)
% Return the interval [i1, i2] such that obj.times(i1) <= t <= obj.times(i2)

if t <= obj.times(1)
  i1 = 1;
  i2 = 1;
  return;
elseif t >= obj.times(end)
  i1 = obj.nTimes;
  i2 = obj.nTimes;
  return;
end

i1 = 1;
i2 = obj.nTimes;

while (i2 - i1 > 1)
  pos = floor((i1 + i2)/2);
  if obj.times(pos) > t
    i2 = pos;
  else
    i1 = pos;
  end
end
end
