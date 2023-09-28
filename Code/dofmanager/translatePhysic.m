function newStr = translatePhysic(str)
%Translate physics string in order to correctly query modelType class
%   Detailed explanation goes here
switch str
    case "Poromechanics"
        newStr = 'Poro';
    case "SPFlow"
        newStr = 'Flow';
    case "Poro"
        newStr = 'Poromechanics';
    case "Flow"
        newStr = 'SPFlow';
end
end

