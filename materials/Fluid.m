classdef Fluid < handle
  % FLUID general fluid class  
    
  properties (Access = private)
    %General properties:
    gamma = [];             %Fluid specific weight 
    beta = [];              %Fluid compressibility 
  end 
    
  methods (Access = public)
    % Class constructor method
    function obj = Fluid(inputString)
      % Calling the function setMaterialParameters to set material
      % parameters
      obj.setMaterialParameters(inputString);
    end
    
    % Function to get material weight
    function gamma = getWeight(obj)
      gamma = obj.gamma;
    end

    % Function to get material compressibility
    function beta = getCompressibility(obj)
      beta = obj.beta;
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
      obj.gamma = params(1);
      obj.beta = params(2);
    end
  end

end
