function str = translatePhysic(str,model)
% Translate physics string in order to correctly query different classes
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
   case "Poro"
      str = "Poromechanics";
end
end