classdef TransvElastic < handle
  % ELASTIC TRANSVERSE ISOTROPIC material class

  properties (Access = private)
    % Elastic modulus in the symmetry plane (x,y)
    E_p = [];
    % Elastic modulus in perpendicular direction (z) 
    E_z = [];
    % Poisson's ratio in the symmetry plane (x,y)
    nu_p = [];
    % Poisson's ratio in perpendicular direction (z) 
    nu_z = [];
  end

  methods (Access = public)
    % Class constructor method
    function obj = TransvElastic(inputString)
      % Calling the function to set object properties 
      obj.setMaterialParameters(inputString);
    end

    % Material stiffness matrix calculation using the object properties
    function D = getStiffnessMatrix(obj, varargin)
      % Elastic moduli ratio
      lambda = obj.E_p/obj.E_z;
      % Constituent matrix
      D = zeros(6,6);
      D(1,1) = obj.E_p*(lambda-(obj.nu_z)^2)/((lambda-lambda*obj.nu_p-2*(obj.nu_z)^2)*(1+obj.nu_p));
      D(1,2) = obj.E_p*(lambda*obj.nu_p+(obj.nu_z)^2)/((lambda-lambda*obj.nu_p-2*(obj.nu_z)^2)*(1+obj.nu_p));
      D(1,3) = obj.E_p*obj.nu_z/(lambda-lambda*obj.nu_p-2*(obj.nu_z)^2);
      D(2,1) = obj.E_p*(lambda*obj.nu_p+(obj.nu_z)^2)/((lambda-lambda*obj.nu_p-2*(obj.nu_z)^2)*(1+obj.nu_p));
      D(2,2) = obj.E_p*(lambda-(obj.nu_z)^2)/((lambda-lambda*obj.nu_p-2*(obj.nu_z)^2)*(1+obj.nu_p));
      D(2,3) = obj.E_p*obj.nu_z/(lambda-lambda*obj.nu_p-2*(obj.nu_z)^2);
      D(3,1) = obj.E_p*obj.nu_z/(lambda-lambda*obj.nu_p-2*(obj.nu_z)^2);
      D(3,2) = obj.E_p*obj.nu_z/(lambda-lambda*obj.nu_p-2*(obj.nu_z)^2);
      D(3,3) = ((1-obj.nu_p)*obj.E_p)/(lambda-lambda*obj.nu_p-2*(obj.nu_z)^2);
      D(4,4) = obj.E_p/(2*(1+obj.nu_z));
      D(5,5) = obj.E_p/(2*(1+obj.nu_z));
      D(6,6) = 0.5*(D(1,1)-D(1,2));
      D = abs(D);
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
      % Object properties are assigned with the same order used for material parameters inside the input file
      obj.E_p = params(1);
      obj.E_z = params(2);
      obj.nu_p = params(3);
      obj.nu_z = params(4);
    end
  end

end
