classdef PorousRock < handle
  %MATERIAL General material class

  properties (Access = private)
    %General properties:
    kx = [];             %Elastic modulus 
    ky = [];            %Poisson ratio 
    kz = [];
    poro = [];
  end

  methods (Access = public)
    function obj = PorousRock(inputString)
      obj.setMaterialParameters(inputString);
    end

    function poro = getPorosity(obj)
      poro = obj.poro;
    end

    function [kx, ky, kz] = getPermeability(obj)
      kx = obj.kx;
      ky = obj.ky;
      kz = obj.kz;
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
      obj.kx = params(1);
      obj.ky = params(2);
      obj.kz = params(3);
      obj.poro = params(4);
    end
  end

end
