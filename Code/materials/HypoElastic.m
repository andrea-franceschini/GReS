classdef HypoElastic < handle
  % HYPOELASTIC ISOTROPIC material class

  properties (Access = private)
    % Poisson's ratio
    nu
    % Coefficients a, b for virgin compressibility 
    a
    b
    % Coefficients a1, b1 for unload/reload compressibility
    a1
    b1
    % Preconsolidation vertical stress
    szmin
    D1
    % M factor (sigmax/sigmaz)
    M
  end

  methods (Access = public)
    % Class constructor method
    function obj = HypoElastic(inputString)
      % Calling the function to set the object properties
      obj.setMaterialParameters(inputString);
      obj.computeConstPart();
    end

    % Material stiffness matrix calculation using the object properties
%     function D = getStiffnessMatrix(obj, varargin)
%       if (nargin ~= 2) 
%         error('Error in calling the HypoElastic/getStiffnessMatrix method - INPUT: (sz)');
%       end
%       % Stiffness matrix
%       % varargin{1} is sz, i.e., the vertical stress 
%       cm = getCompressibility(obj, varargin{1});
%       D = (1/cm)*obj.D1;
%     end
    function D = getStiffnessMatrix(obj, sz)
%       if (nargin ~= 2) 
%         error('Error in calling the HypoElastic/getStiffnessMatrix method - INPUT: (sz)');
%       end
      % Stiffness matrix
      % varargin{1} is sz, i.e., the vertical stress 
      cM = getRockCompressibility(obj, sz);
      D = obj.D1.*reshape(1./cM,1,1,[]);
    end
    
    % Method that returns the M factor
    function m = getMFactor(obj)
      m = obj.M;
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
      if nEntry ~= 6
        error('Wrong number of input parameters in material %s',block(1));
      end
      params = str2double(strParams);
      % Assign object properties
      obj.nu = params(1);
      obj.a = params(2);
      obj.b = params(3);
      obj.a1 = params(4);
      obj.b1 = params(5);
      obj.szmin = params(6);
      %
      % Compute the M factor
      obj.M = obj.nu/(1-obj.nu);
    end
    
    function computeConstPart(obj)
      % Stiffness matrix
      obj.D1 = zeros(6);
      obj.D1([1 8 15]) = 1;
      obj.D1([2 3 7 9 13 14]) = obj.nu/(1-obj.nu);
      obj.D1([22 29 36]) = (1-2*obj.nu)/(2*(1-obj.nu));
    end

    % Compressibility calculation
    function cM = getRockCompressibility(obj, sz)
      if sz <= obj.szmin
        % Loading path
        cM = (obj.a).*(abs(sz)).^(obj.b);
      else
        % Unloading/reloading path
        cM = (obj.a1).*(abs(sz)).^(obj.b1);
      end
    end
  end
end