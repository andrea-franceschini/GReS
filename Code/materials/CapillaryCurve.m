classdef CapillaryCurve < handle
  % ELASTIC ISOTROPIC material class

  properties (Access = private)
    tab
  end

  methods (Access = public)
    % Class constructor method
    function obj = CapillaryCurve(fID)
      % Calling the function to set the object properties 
      obj.readCapillaryCurve(fID);
    end
  end

  methods (Access = private)
    function readCapillaryCurve(obj,fID)
      obj.tab = load(fID);
    end
  end
end