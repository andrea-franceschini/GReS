classdef structGrid
  %UNTITLED3 Summary of this class goes here
  %   Detailed explanation goes here

  properties
    arrayX
    arrayY
    arrayZ
    divisions (1,3)
  end

  methods
    function obj = structGrid(data)
      %structGrid Construct an instance of this class
      %   Detailed explanation goes here
      switch lower(data.type)
        case "classic"
          obj.constructorClassic(data);
        otherwise
      end
    end

    function outputArg = method1(obj,inputArg)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      outputArg = obj.Property1 + inputArg;
    end
  end

  methods (Access = private)
    function obj = constructorClassic(obj,data)
      obj.divisions = str2num(data.division);
      dim = str2num(data.size);
      obj.arrayX = linspace(0,dim(1),obj.divisions(1));
      obj.arrayY = linspace(0,dim(2),obj.divisions(2));
      obj.arrayZ = linspace(0,dim(3),obj.divisions(3));
    end

    
  end

end