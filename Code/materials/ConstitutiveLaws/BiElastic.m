classdef BiElastic < Elastic
  % A bielastic constitutive law class where different compressibility is
  % assumed during virgin loading or unloading/reloading

  properties (Access = public)
    % scaling of cm between loading and unloading 
    r
    % preconsolidation stress
    q
    % loading state (0 is virigin, 1 is unloading/realoding)
    loadState

  end

  methods (Access = public)
    % Class constructor method
    function obj = BiElastic(varargin)
          readMaterialParameters(obj,varargin{:});
    end
    
    function [status] = initializeStatus(obj, sigma)
      nptGauss = size(sigma,1);
      status = zeros(nptGauss,6);
    end
    %
    % Material stiffness matrix calculation using the object properties
    function [DAll, sigmaOut, status] = getStiffnessMatrix(obj, sigmaIn, epsilon, dt, status, cellID)
      
      nptGauss = size(sigmaIn,1);
      D = getElasticTensor(obj,cellID);

      % we assume compression
      p = abs(sum(sigma(1:3))/3);

      obj.loadState(cellID) = 0;

      if p < obj.q(cellID)
        obj.loadState(cellID) = 1;
        D = obj.r * D;
      end

      sigmaOut = sigmaIn + epsilon*D;

      % neglect piecewise transition from reloading to virgin loading

      DAll = repmat(D,[1, 1, nptGauss]);
    end
    %

 
    % Method that returns the M factor
    function m = getMFactor(obj)
      m = obj.M;
    end
    
    % Get vertical compressibility
    function cM = getRockCompressibility(obj,cellIds)

      % scale factor for rock compressibility
      scale = 1 - (1/obj.r) .* (obj.loadState(cellIds) == 1);  
      cM = scale * obj.cM;

    end
  end
   
  methods (Access = private)
    % Assigning material parameters (check also the Materials class)
    % to object properties
    function readMaterialParameters(obj,varargin)

      default = struct("youngModulus",[],...
                       "poissonRatio",[],...
                       "");

      params = readInput(default,varargin{:});
      
       %
       obj.E = params.youngModulus;
       obj.nu = params.poissonRatio;
       %
       % Compute the M factor
       obj.M = obj.nu/(1-obj.nu);
       %
       % Compute vertical compressibility
       obj.cM =  (1+obj.nu).*(1-2*obj.nu)./(obj.E.*(1-obj.nu));

    end
 
    % function readTabMaterialParameters(obj,fID,fileName,mesh)
    %    % young modulus
    %    youngModFile = readToken(fID,fileName);
    %    poissonRatioFile = readToken(fID,fileName);
    %    obj.E = setTabularParams(youngModFile,mesh);
    %    obj.nu = setTabularParams(poissonRatioFile,mesh);
    %    % Compute the M factor
    %    obj.M = obj.nu./(1-obj.nu);
    %    %
    %    % Compute vertical compressibility
    %    obj.cM = (1+obj.nu).*(1-2*obj.nu)./(obj.E.*(1-obj.nu));
    % end
  end
end