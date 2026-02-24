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
    % flag for tabular input format
    isTabular = false
  end

  methods (Access = public)
    % Class constructor method
    function obj = Elastic(varargin)
          readMaterialParameters(obj,varargin{:});
    end
    
    function [status] = initializeStatus(obj, sigma)
      nptGauss = size(sigma,1);
      status = zeros(nptGauss,2);
    end
    %
    % Material stiffness matrix calculation using the object properties
    function [DAll, sigmaOut, status] = getStiffnessMatrix(obj, sigmaIn, epsilon, dt, status, cellID)
      nptGauss = size(sigmaIn,1);
      D = getElasticTensor(obj,cellID);
      sigmaOut = sigmaIn + epsilon*D;
      DAll = repmat(D,[1, 1, nptGauss]);
    end
    %

    function D = getElasticTensor(obj,cID)
      % elastic tensor in engineering Voigt notation
       D = zeros(6);
       if obj.isTabular
          pois = obj.nu(cID);
          D([1 8 15]) = 1-pois;
          D([2 3 7 9 13 14]) = pois;
          D([22 29 36]) = (1-2*pois)/2;
          D = obj.E(cID)/((1+pois)*(1-2*pois))*D;
       else
          D([1 8 15]) = 1-obj.nu;
          D([2 3 7 9 13 14]) = obj.nu;
          D([22 29 36]) = (1-2*obj.nu)/2;
          D = obj.E/((1+obj.nu)*(1-2*obj.nu))*D;
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
    function readMaterialParameters(obj,varargin)
      
      if isfield(inputStruct,"Elastic")
        inputStruct = inputStruct.Elastic;
      end
       %
       obj.E = getXMLData(inputStruct,[],"youngModulus");
       obj.nu = getXMLData(inputStruct,[],"poissonRatio");
       %
       % Compute the M factor
       obj.M = obj.nu/(1-obj.nu);
       %
       % Compute vertical compressibility
       obj.cM = 0.;
       % Compute constitutive matrix
    end
 
    function readTabMaterialParameters(obj,fID,fileName,mesh)
       % young modulus
       youngModFile = readToken(fID,fileName);
       poissonRatioFile = readToken(fID,fileName);
       obj.E = setTabularParams(youngModFile,mesh);
       obj.nu = setTabularParams(poissonRatioFile,mesh);
       % Compute the M factor
       obj.M = obj.nu./(1-obj.nu);
       %
       % Compute vertical compressibility
       obj.cM = (1+obj.nu).*(1-2*obj.nu)./(obj.E.*(1-obj.nu));
    end
  end
end