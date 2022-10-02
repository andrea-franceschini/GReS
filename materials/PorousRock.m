classdef PorousRock < handle
  % POROUS ROCK material class

  properties (Access = private)
    %General properties:
    kx = [];             %Permeability in x
    ky = [];             %Permeability in y 
    kz = [];             %Permeability in z
    poro = [];           %Porosity
  end

  methods (Access = public)
      % Class constructor method
    function obj = PorousRock(inputString)
      % Calling the function setMaterialParameters to set material
      % parameters
      obj.setMaterialParameters(inputString);
    end

    % Function to get material porosity
    function poro = getPorosity(obj)
      poro = obj.poro;
    end

    % Function to get material permeability
    function [kx, ky, kz] = getPermeability(obj)
      kx = obj.kx;
      ky = obj.ky;
      kz = obj.kz;
    end
  end

  methods (Access = private)
    % Function that assigns the material data the object properties
    % coming the cell "data" created inside the class Materials
    function setMaterialParameters(obj, inputString)
      words = strsplit(inputString, ' ');
      params = zeros(length(words),1);
      k = 0;
      for i = 1 : length(words)
        if (length(words{i}) > 0)
          k = k + 1;
          params(k) = sscanf(words{i}, '%e');
        end
      end
      % Object properties are assigned with the same order used for the 
      % material parameters inside the input file
      obj.kx = params(1);
      obj.ky = params(2);
      obj.kz = params(3);
      obj.poro = params(4);
    end
  end

end
