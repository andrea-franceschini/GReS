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
    function obj = Fluid(inputString)
      % Calling the function to set the material parameters
      obj.setMaterialParameters(inputString);
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
    function setMaterialParameters(obj, block)
      % Preliminary check on the number of rows in each material block
      % and the number of parameters
      nEntry = size(block,1);
      if nEntry ~= 2
        error('Wrong number of input rows in material %s', block(1));
      end
      strParams = strsplit(block(2));
      nEntry = size(strParams,2);
      if nEntry ~= 3
        error('Wrong number of input parameters in material %s',block(1));
      end
      params = str2double(strParams);
      % Assign object properties
      obj.gamma = params(1);
      obj.beta = params(2);
      obj.mu = params(3);
    end
  end

end