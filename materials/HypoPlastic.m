classdef HypoPlastic < handle
  % Hypoplastic material class

  properties (Access = private)
    % Poisson's ratio
    nu = [];
    % a, b, in cm = a * sz^b
    a = [];
    b = [];
    % Initial void index
    e0 = [];
    % Compression coefficient
    Cc = [];
    % Recompression index
    Cr = [];
  end

  methods (Access = public)
    function obj = HypoPlastic(inputString)
      obj.setMaterialParameters(inputString);
    end

    function D = getStiffnessMatrix(obj, varargin)
      if (nargin < 2)
        error('Missing sz in HypoPlastic/getStiffnessMatrix');
      end
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
      D = 1/cm*((1+obj.nu)*(1-2*obj.nu))*D;
    end
  end

  methods (Access = private)
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
      obj.nu = params(1);
      obj.a = params(2);
      obj.b = params(3);
      obj.e0 = params(4);
      obj.Cc = params(5);
      obj.Cr = params(6);
    end

    function cm = getCompressibility(obj, sz)
      cm = obj.a * sz.^obj.b;
    end
  end
end
