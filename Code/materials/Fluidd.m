classdef Fluidd < MatEntry
  % FLUID general fluid class

  properties (Access = private)
    %General properties:
    name
    gamma;             % Fluid specific weight
    beta;              % Fluid compressibility
    mu;                % Fluid dynamic viscosity
  end

  methods (Access = public)
    function obj = Fluidd(varargin)
      if ~isempty(varargin)
        % Calling the function to set the material parameters
        obj.readMaterialParameters(varargin{:});
      end
    end

    function gamma = getSpecificWeight(obj)
      %GETFLUIDSPECWIEGHT Function to get fluid specific gravity
      gamma = obj.gamma;
    end

    function beta = getFluidCompressibility(obj)
      % GETFLUIDCOMPRESSIBILITY Function to get fluid compressibility
      beta = obj.beta;
    end

    function mu = getDynViscosity(obj)
      %GETDYNVISCOSITY Function to get fluid dynamic viscosity
      mu = obj.mu;
    end

    function setName(obj,name)
      obj.name = name;
    end

    function materialUpdate(obj,prop,varargin)
    end

    % compute the jacobian and the rhs
    function getProperty(obj,prop,varargin)
    end

  end

  methods (Access = private)
    function readMaterialParameters(obj,varargin)

      default = struct('specificWeight',0.0,...
        'compressibility',0.0,...
        'dynamicViscosity',1e-6); % viscosity of water in kPa*s

      params = readInput(default,varargin{:});

      obj.registerProp('specificWeight',params.specificWeight);
      obj.registerProp('compressibility',params.compressibility);
      obj.registerProp('dynamicViscosity',params.dynamicViscosity);
    end
  end

  methods (Static)
    function out = getProps()
      out = cell("SpecificWeight","Compressibility","DynamicViscosity");
    end
    
  end
end