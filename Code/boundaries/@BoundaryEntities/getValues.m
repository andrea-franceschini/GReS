function vals = getValues(obj, t)

if (obj.nTimes==0)
  % value must a function handle depending on time
  tScale = str2func(['@(t)', char(obj.bcData.time)]);
  if obj.readSetFlag
    obj.availVals(:,1) = readDataSet(obj);
    obj.readSetFlag = false;
  end
  vals = tScale(t)*obj.availVals(:,1); % scale by the time function
elseif (obj.nTimes == 1)
  if obj.readSetFlag
    obj.availVals(:,1) = readDataSet(obj);
    obj.readSetFlag = false;
  end
  vals = obj.availVals(:,1);
else
  [i1, i2] = bin_search(obj, t);
  if i1 == i2 % edge case
    if all([i1 i2] == 1)
      obj.availSteps(1) = i1;
      obj.availVals(:,1) = readDataSet(obj,1);
      vals = obj.availVals(:,1);
    elseif all([i1 i2] == obj.nTimes)
      obj.availVals(:,2) = readDataSet(obj,obj.nTimes);
      vals = obj.availVals(:,2);
    end
  else
    p1 = find(obj.availSteps == i1);
    p2 = find(obj.availSteps == i2);
    if (isempty(p1) && isempty(p2))
      p1 = 1;
      obj.availVals(:,1) = readDataSet(obj,i1);
      obj.availSteps(1) = i1;
      p2 = 2;
      obj.availVals(:,2) = readDataSet(obj,i2);
      obj.availSteps(2) = i2;
    elseif (~isempty(p1) && isempty(p2))
      p2 = 3 - p1;
      obj.availVals(:,p2) = readDataSet(obj,i2);
      obj.availSteps(p2) = i2;
    elseif (~isempty(p2) && isempty(p1))
      p1 = 3 - p2;
      obj.availVals(:,p1) = readDataSet(obj,i1);
      obj.availSteps(p1) = i1;
    end
    fac = (t - obj.times(i1)) / (obj.times(i2) - obj.times(i1));
    vals = fac*(obj.availVals(:,p2) - obj.availVals(:,p1)) + obj.availVals(:,p1);
  end

end

%vals = vals(obj.isActiveEntity);
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
