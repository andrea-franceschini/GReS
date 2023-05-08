classdef Elastic < handle
  % ELASTIC ISOTROPIC material class

  properties (Access = public)
    % Elastic modulus
    E
    % Poisson's ratio
    nu
    % M factor (sigmax/sigmaz)
    M
    % Vertical compressibility Cm
    cM
  end

  methods (Access = public)
    % Class constructor method
    function obj = Elastic(inputString)
      % Calling the function to set the object properties 
      obj.setMaterialParameters(inputString);
    end
    %
    % Material stiffness matrix calculation using the object properties
    function D = getStiffnessMatrix(obj)
      % Stiffness matrix
      D = zeros(6);
      D([1 8 15]) = 1-obj.nu;
      D([2 3 7 9 13 14]) = obj.nu;
      D([22 29 36]) = (1-2*obj.nu)/2;
      D = obj.E/((1+obj.nu)*(1-2*obj.nu))*D;
    end
    
    % Method that returns the M factor
    function m = getMFactor(obj)
      m = obj.M;
    end
    
    % Get vertical compressibility
    function cM = getRockCompressibility(obj)
      cM = obj.cM;
    end
  end

  methods (Access = private)
    % Assigning material parameters (check also the Materials class)
    % to object properties
    function setMaterialParameters(obj,block)
      % Preliminary check on the number of rows in each material block
      % and the number of parameters
      nEntry = size(block,1);
      if nEntry ~= 2
        error('Wrong number of input rows in material %s',block(1));
      end
      strParams = strsplit(block(2));
      nEntry = size(strParams,2);
      if nEntry ~= 2
        error('Wrong number of input parameters in material %s',block(1));
      end
      params = str2double(strParams);
      % Assign object properties
      obj.E = params(1);
      obj.nu = params(2);
      %
      % Compute the M factor
      obj.M = obj.nu/(1-obj.nu);
      %
      % Compute vertical compressibility
      obj.cM = (1+obj.nu)*(1-2*obj.nu)/(obj.E*(1-obj.nu));
    end
  end
end