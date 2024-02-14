function str = translatePhysic(str,model)
%Translate physics string in order to correctly query different classes
%   Detailed explanation goes here
switch str
    case "Poromechanics"
        str = "Poro";
    case "SPFlow"
        str = "Flow";
    case "VSFlow"
        str = "Flow";
    case "Flow"
        if isSinglePhaseFlow(model)
            str = "SPFlow";
        elseif isVariabSatFlow(model)
            str = "VSFlow";
        end
end
end