classdef MaterialPool
  %MATERIALPOOL Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    data = struct()
  end
  
  methods
    function obj = MaterialPool(inputArg1,inputArg2)
      %MATERIALPOOL Construct an instance of this class
      %   Detailed explanation goes here
      obj.Property1 = inputArg1 + inputArg2;
    end
    
    function outputArg = getProperty(obj,propName,varargin)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here

      % varargin: cell Id
      prop
    end
  end
end

