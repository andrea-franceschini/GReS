classdef UserDefElastic < handle
  % ELASTIC ISOTROPIC material class

  properties (Access = public)
    % UserDefElastic stiffness array
    stiffnesses
  end

  methods (Access = public)
    % Class constructor method
    function obj = UserDefElastic(varargin)
      % Unused: just for compatibility!
    end

    % setter
    function obj = setStiffnesses(obj, stiffnesses)
      obj.stiffnesses = stiffnesses;
    end
    
    function [status] = initializeStatus(obj, sigma)
      nptGauss = size(sigma,1);
      status = zeros(nptGauss,2);
    end
    %
    % Material stiffness matrix calculation using the object properties
    function [DAll, sigmaOut, status] = getStiffnessMatrix(obj, sigmaIn, epsilon, dt, status, cellID)
      nptGauss = size(sigmaIn,1);
      D = obj.stiffnesses{cellID};
      sigmaOut = sigmaIn + epsilon*D;
      DAll = repmat(D,[1, 1, nptGauss]);
    end
    %
  end
end
