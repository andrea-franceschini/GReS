function sol = findzeros(f,range,err)
if nargin < 3
  err = 1e-3;
end
sol = vpasolve(f,range);
if(isempty(sol))
  return
else
  lowLimit = sol-err;
  highLimit = sol+err;
  temp = findzeros(f,[range(1) lowLimit],1);
  if ~isempty(temp)
    sol = sort([sol temp]);
  end
  temp = findzeros(f,[highLimit range(2)],1);
  if ~isempty(temp)
    sol = sort([sol temp]);
  end
  return
end
end