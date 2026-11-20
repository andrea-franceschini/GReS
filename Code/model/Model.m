classdef Model < handle
  % GReS model manager
  
  properties

    solutionStrategy
    domains
    interfaces

    J       
    res

    simParam
    output


  end
  
  methods
    function obj = Model(inputArg1,inputArg2)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      obj.Property1 = inputArg1 + inputArg2;
    end
    
    function outputArg = method1(obj,inputArg)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      outputArg = obj.Property1 + inputArg;
    end
  end
end

