classdef Elastic < handle
  % ELASTIC ISOTROPIC material class

  properties (Access = private)
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
    function obj = Elastic(fID, matFileName)
      % Calling the function to set the object properties 
      obj.readMaterialParameters(fID, matFileName);
    end

    function [status] = initializeStatus(obj, sigma)
      nptGauss = size(sigma,1);
      status = zeros(nptGauss,2);
    end
    %
    % Material stiffness matrix calculation using the object properties
    function [DAll, sigmaOut, status] = getStiffnessMatrix(obj, sigmaIn, eps, dt, status)
      nptGauss = size(sigmaIn,1);
      sigmaOut = zeros(nptGauss,6);
      % Stiffness matrix
      DAll = zeros(6,6,nptGauss);
      D = zeros(6,6);
      D([1 8 15]) = 1-obj.nu;
      D([2 3 7 9 13 14]) = obj.nu;
      D([22 29 36]) = (1-2*obj.nu)/2;
      D = obj.E/((1+obj.nu)*(1-2*obj.nu))*D;
      for i = 1 : nptGauss
        sigmaOut(i,:) = sigmaIn(i,:) + eps(i,:)*D;
        DAll(:,:,i) = D;
      end
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
    function readMaterialParameters(obj, fID, matFileName)
      tmpVec = readDataInLine(fID, matFileName, 2);
      %
      obj.E = tmpVec(1);
      obj.nu = tmpVec(2);
      %
      % Compute the M factor
      obj.M = obj.nu/(1-obj.nu);
      %
      % Compute vertical compressibility
      obj.cM = (1+obj.nu)*(1-2*obj.nu)/(obj.E*(1-obj.nu));
    end
  end
end
