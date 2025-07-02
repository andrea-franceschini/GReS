classdef Fluid < handle
  % FLUID general fluid class  

  properties (Access = private)
    %General properties:
    gamma;             % Fluid specific weight 
    beta;              % Fluid compressibility
    mu;                % Fluid dynamic viscosity
  end 

  methods (Access = public)
    % Class constructor method
    function obj = Fluid(fID, matFileName)
      % Calling the function to set the material parameters
      obj.readMaterialParameters(fID, matFileName);
    end

    % Function to get fluid specific gravity
    function gamma = getFluidSpecWeight(obj)
      gamma = obj.gamma;
    end

    % Function to get fluid compressibility
    function beta = getFluidCompressibility(obj)
      beta = obj.beta;
    end
    
    % Function to get fluid dynamic viscosity
    function mu = getDynViscosity(obj)
      mu = obj.mu;
    end
  end

  methods (Access = private)
    % Assigning material parameters (check also the Materials class)
    % to object properties
    function readMaterialParameters(obj, fID, matFileName)
      tmpVec = readDataInLine(fID, matFileName, 3);
      % Assign object properties
      obj.gamma = tmpVec(1);
      obj.beta = tmpVec(2);
      obj.mu = tmpVec(3);
    end
  end

end