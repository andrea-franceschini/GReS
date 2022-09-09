classdef HypoElastic < handle
  % Hypoelastic isotropic material class

  properties (Access = private)
    % Poisson's ratio
    nu = [];
    % a, b, in cm = a * sz^b
    a = [];
    b = [];
    % a1, b1, in cm = a1 * sz^b1
    a1 = [];
    b1 = [];
    % Lowest vertical stress
    szmin = [];
  end

  methods (Access = public)
      % Class constructor method
    function obj = HypoElastic(inputString)
      % Calling the function setMaterialParameters to set material
      % parameters
      obj.setMaterialParameters(inputString);
    end

    % Function calculating the stiffness matrix using the class properties
    function D = getStiffnessMatrix(obj, varargin)
      % Setting the minimum number of input values in order to calcolate the
      % constituent matrix
      if (nargin < 2) 
        error('Missing sz in HypoPlastic/getStiffnessMatrix');
      end
      % Setting z-stress as the first getStiffnessMatrix input value
      sz = varargin{1};
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
      cm = getCompressibility(obj, sz);
      D = (1/cm)*D;
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
      obj.nu = params(1);
      obj.a = params(2);
      obj.b = params(3);
      obj.a1 = params(4);
      obj.b1 = params(5);
      obj.szmin = params(6);
    end

    % Function for compressibility calculation
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
