classdef ModCamClay < handle
  % Modified Cam Clay material class

  properties (Access = private)
    %Stiffness parameters:
    lambda = [];             %Isotropic compression index
    k = [];                  %Isotropic swelling index 
    nu_ur = [];              %Unloading-reloading Poisson ratio
    e0 = [];                 %Initial void ratio
    %Strenght parameters:
    M = [];                  %Tangent of the critical state line
    K0 = [];                 %Coefficient of lateral stress in normal consolidation
  end

  methods (Access = public)
    function obj = ModCamClay(inputString)
      obj.setMaterialParameters(inputString);
    end

%     function lambda = getStiffness(obj)
%       poro = obj.poro;
%     end
% 
%     function [kx, ky, kz] = getPermeability(obj)
%       kx = obj.kx;
%       ky = obj.ky;
%       kz = obj.kz;
%     end
%   end

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
%       obj.lambda = params(1);
%       obj.k = params(2);
%       obj.nu_ur = params(3);
%       obj.e0 = params(4);
    end
  end

end
