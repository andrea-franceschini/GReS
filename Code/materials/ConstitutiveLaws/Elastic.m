classdef Elastic < ConstitutiveLaw
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
    function obj = Elastic(varargin)
          readMaterialParameters(obj,varargin{:});
    end
    
    function initializeStatus(obj, state, gaussId) %#ok<INUSD>

      % no state variables in Linear elastic law


    end


    %
    % Material stiffness matrix calculation using the object properties
    function [sigmaOut,DAll] = constitutiveUpdate(obj, cellId, sigmaIn, epsilon, varargin)

      nptGauss = obj.loc2gp(cellId(1),2);

      D = getElasticTensor(obj);
      sigmaOut = sigmaIn + epsilon*D;

      % all cells have the same elastic properties
      DAll = repmat(D,[1, 1, nptGauss,numel(cellId)]);

    end
    %

    function D = getElasticTensor(obj)
      % elastic tensor in engineering Voigt notation
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
    function cM = getRockCompressibility(obj,cId) %#ok<INUSD>
      cM = obj.cM;
    end
  end
   
  methods (Access = protected)

    % Assigning material parameters (check also the Materials class)
    % to object properties
    function readMaterialParameters(obj,varargin)

      default = struct("youngModulus",[],...
                       "poissonRatio",[]);

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

  end


  methods (Static)

    function out = isLinear()

      out = true;

    end
  end
end