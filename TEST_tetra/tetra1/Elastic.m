classdef Elastic < handle
  % ELASTIC ISOTROPIC material class

  properties (Access = private)
    % Elastic modulus
    E = [];
    % Poisson's ratio
    nu = [];
  end

  methods (Access = public)
    % Class constructor method
    function obj = Elastic(inputString)
      % Calling the function to set object properties 
      obj.setMaterialParameters(inputString);
    end

    % Material stiffness matrix calculation using the object properties
    function D = getStiffnessMatrix(obj, varargin)
      % Constituent matrix
      D = zeros(6,6);
      D(1,1) = 1-obj.nu;
      D(1,2) = obj.nu;
      D(1,3) = obj.nu;
      D(2,1) = obj.nu;
      D(2,2) = 1-obj.nu;
      D(2,3) = obj.nu;
      D(3,1) = obj.nu;
      D(3,2) = obj.nu;
      D(3,3) = 1-obj.nu;
      D(4,4) = (1-2*obj.nu)/2;
      D(5,5) = (1-2*obj.nu)/2;
      D(6,6) = (1-2*obj.nu)/2;
      D = obj.E/((1+obj.nu)*(1-2*obj.nu))*D;
    end
  end

  methods (Access = private)
    % Assigning material parameters (read inside the class Materials) to object properties
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
      % Object properties are assigned with the same order used for
      % material parameters inside the input file
      obj.E = params(1);
      obj.nu = params(2);
    end
  end

end