function topol = build_topol(nodes)
%BUILD_TOPOL Summary of this function goes here
%   Detailed explanation goes here
topol = 1:length(nodes);
topol = [topol(1) repelem(topol(2:end-1),2), topol(end)];
topol = (reshape(topol, 2, []))';
end