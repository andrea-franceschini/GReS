function str = translatePhysic(str)
%Translate physics string in order to correctly query different classes
%   Detailed explanation goes here
if strcmp(str, "Poromechanics")
    str = "Poro";
end
if strcmp(str, "SPFlow")
    str = "Flow";
end
end

