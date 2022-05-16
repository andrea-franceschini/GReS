classdef TransvElastic < handle
  % Elastic anisotropic (transversely isotropic) material class

  properties (Access = private)
    % Elastic modulus in the x-y symmetry plane (plane1)
    E_p = [];
    % Elastic modulus in z-direction (perpendicular to plane1) 
    E_z = [];
    % Poisson's ratio in the x-y symmetry plane (plane1)
    nu_p = [];
    % Poisson's ratio in z-direction (perpendicular to plane1) 
    nu_z = [];
  end

  methods (Access = public)
      % Class constructor method
    function obj = TransvElastic(inputString)
      % Calling the function setMaterialParameters to set material
      % parameters
      obj.setMaterialParameters(inputString);
    end

    % Function calculating the stiffness matrix using the class properties
    function D = getStiffnessMatrix(obj, varargin)
      delta = (1+obj.nu_p)*(1-obj.nu_p-2*(obj.nu_z)^2)/((obj.E_p)^2+obj.E_z);
      % Constituent matrix
      D = zeros(6,6);
      D(1,1) = (1-(obj.nu_z)^2)/(obj.E_p*obj.E_z*delta);
      D(1,2) = (obj.nu_p+(obj.nu_z)^2)/(obj.E_p*obj.E_z*delta);
      D(1,3) = (obj.nu_z+obj.nu_z*obj.nu_p)/(obj.E_p*obj.E_z*delta);
      D(2,1) = (obj.nu_p+(obj.nu_z)^2)/(obj.E_p*obj.E_z*delta);
      D(2,2) = (1-(obj.nu_z)^2)/(obj.E_p*obj.E_z*delta);
      D(2,3) = (obj.nu_z+obj.nu_z*obj.nu_p)/(obj.E_p*obj.E_z*delta);
      D(3,1) = (obj.nu_z+obj.nu_z*obj.nu_p)/(obj.E_p*obj.E_z*delta);
      D(3,2) = (obj.nu_z+obj.nu_z*obj.nu_p)/(obj.E_p*obj.E_z*delta);
      D(3,3) = 1-(obj.nu_p)^2/((obj.E_p)^2*delta);
      D(4,4) = obj.E_z/(1+obj.nu_z);
      D(5,5) = obj.E_z/(1+obj.nu_z);
      D(6,6) = obj.E_p/(1+obj.nu_p);
    end
  end

  methods (Access = private)
      % Function that set the material parameters coming from "data"
      % (Materials) inside the vector "params"
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
      % Object properties are assigned with the same order used in the input
      % file
      obj.E_p = params(1);
      obj.E_z = params(2);
      obj.nu_p = params(3);
      obj.nu_z = params(4);
    end
  end

end
