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
    function obj = HypoElastic(fID, matFileName)
      % Calling the function to set the object properties
      obj.readMaterialParameters(fID, matFileName);
      obj.computeConstPart();
    end
    
        function [status] = initializeStatus(obj, sigma)
      nptGauss = size(sigma,1);
      status = zeros(nptGauss,2);
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
    function [DAll, sigmaOut, status] = getStiffnessMatrix(obj, sigmaIn, epsilon, dt, status, t)
%       if (nargin ~= 2) 
%         error('Error in calling the HypoElastic/getStiffnessMatrix method - INPUT: (sz)');
%       end
      % Stiffness matrix
      % varargin{1} is sz, i.e., the vertical stress
      nptGauss = size(sigmaIn,1);
      cM = getRockCompressibility(obj, sigmaIn, epsilon, t); %medium value in the element
      D = obj.D1.*reshape(1./cM,1,1,[]);
      sigmaOut = sigmaIn + epsilon*D;
      DAll = repmat(D,[1, 1, nptGauss]);
    end
    
    % Method that returns the M factor
    function m = getMFactor(obj)
      m = obj.M;
    end
  end

  methods (Access = private)
    % Assigning material parameters (check also the Materials class)
    % to object properties
    function readMaterialParameters(obj, fID, matFileName)
      tmpVec = readDataInLine(fID, matFileName, 6);
      % Assign object properties
      obj.nu = tmpVec(1);
      obj.a = tmpVec(2);
      obj.b = tmpVec(3);
      obj.a1 = tmpVec(4);
      obj.b1 = tmpVec(5);
      obj.szmin = tmpVec(6);
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
    function cM = getRockCompressibility(obj, sigmaIn, epsilon, t)
        sz = mean(sigmaIn(:,3));
      if sz<obj.szmin
        % Loading path
        cM = (obj.a).*(abs(sz)).^(obj.b);
      else
        % Unloading/reloading path
        cM = (obj.a1).*(abs(sz)).^(obj.b1);
      end
    end
  end
end