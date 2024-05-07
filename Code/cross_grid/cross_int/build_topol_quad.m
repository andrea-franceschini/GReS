function topol = build_topol_quad(nodes)
%BUILD_TOPOL Summary of this function goes here
%   Detailed explanation goes here
nElems = 0.5*(length(nodes)-1);
assert(mod(nElems,1)==0,['Number of nodes not applicable to quadratic elements: ' ...
    'input must be an odd number > 1']);
topol = zeros(nElems,3);
k = 1;
for i = 1:nElems
    topol(i,:) = [k k+1 k+2];
    k = k+2;
end
end

