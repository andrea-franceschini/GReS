classdef HypoElastic < handle
  % HYPOELASTIC ISOTROPIC material class

  properties (Access = private)
    % Poisson's ratio
    nu = [];
    % Coefficients a, b for virgin compressibility 
    a = [];
    b = [];
    % Coefficients a1, b1 for unload/reload compressibility
    a1 = [];
    b1 = [];
    % Preconsolidation vertical stress
    szmin = [];
  end

  methods (Access = public)
    % Class constructor method
    function obj = HypoElastic(inputString)
      % Calling the function to set object properties
      obj.setMaterialParameters(inputString);
    end

    % Material stiffness matrix calculation using the object properties
    function D = getStiffnessMatrix(obj, varargin)
      if (nargin < 2) 
        error('Missing sz in HypoElastic/getStiffnessMatrix');
      end
      % vertical stress = first 'getStiffnessMatrix' input value
      sz = varargin{1};
      % Constituent matrix
      D = zeros(6,6);  
      D(1,1) = 1;
      D(1,2) = obj.nu/(1-obj.nu);
      D(1,3) = obj.nu/(1-obj.nu);
      D(2,1) = obj.nu/(1-obj.nu);
      D(2,2) = 1;
      D(2,3) = obj.nu/(1-obj.nu);
      D(3,1) = obj.nu/(1-obj.nu);
      D(3,2) = obj.nu/(1-obj.nu);
      D(3,3) = 1;
      D(4,4) = (1-2*obj.nu)/(2*(1-obj.nu));
      D(5,5) = (1-2*obj.nu)/(2*(1-obj.nu));
      D(6,6) = (1-2*obj.nu)/(2*(1-obj.nu));
      cm = getCompressibility(obj, sz);
      D = (1/cm)*D;
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
      obj.nu = params(1);
      obj.a = params(2);
      obj.b = params(3);
      obj.a1 = params(4);
      obj.b1 = params(5);
      obj.szmin = params(6);
    end

    % Compressibility calculation
    function cm = getCompressibility(obj, sz)
       if sz <= obj.szmin
          % Loading path
          cm = (obj.a) * (abs(sz))^(obj.b);
       else
          % Unloading/reloading path
          cm = (obj.a1) * (abs(sz))^(obj.b1);
       end
    end
    
  end
end
