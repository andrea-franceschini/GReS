function physics = rename(physics)
%Rename physics string in order to correctly query different classes
%   Detailed explanation goes here
physics(physics == "Poromechanics") = "Poro";
physics(physics == "SPFlow") = "SPFlow";
physics(physics == "VSFlow") = "VSFlow";
end