function out = compareSet(setA,setB)
% fast checking for current and previous active set
out = true;
if isempty(setA) && isempty(setB)
   out = false;
   return
end
if numel(setA) ~= numel(setB)
   return
elseif any(setA~=setB)
   return
end
out = false;
end

