classdef PorousRock < handle
  % POROUS ROCK material class

  properties (Access = private)
    %General properties:
    KVec                 % Vector of permeabilities
                         % (upper triangular part ordered column-wise)
    poro                 % Porosity
%     alpha                % Rock compressibility (can be replaced by the oedometer test compressibility Cm)
    specGrav             % Specific gravity of rock
  end

  methods (Access = public)
    % Class constructor method
    function obj = PorousRock(inputString)
      % Calling the function to set the object properties
      obj.setMaterialParameters(inputString);
    end

    % Function to get material porosity
    function poro = getPorosity(obj)
      poro = obj.poro;
    end

    % Function to get material permeability
    function K = getPermMatrix(obj)
      if length(obj.KVec) == 3
        K = diag(obj.KVec);
      else
        K = [obj.KVec(1) obj.KVec(2) obj.KVec(4);
             obj.KVec(2) obj.KVec(3) obj.KVec(5);
             obj.KVec(4) obj.KVec(5) obj.KVec(6)];
      end
    end
    
    % Function to get rock compressibility
%     function a = getRockCompressibility(obj)
%       a = obj.alpha;
%     end
  end

  methods (Access = private)
    % Assigning material parameters (check also the Materials class)
    % to object properties
    function setMaterialParameters(obj,block)
      % Preliminary check on the number of rows in each material block
      % and the number of parameters
      nEntry = size(block,1);
      if nEntry ~= 5
        error('Wrong number of input rows in material %s',block(1));
      end
      KTmp = zeros(6,1);
      for i=2:5
        strParams = strsplit(block(i));
        nEntry = size(strParams,2);
        err = false;
        switch i
          case 3
            if nEntry ~= 3; err = true; end
          case {2 4}
            if nEntry ~= 2; err = true; end
          case 5
            if nEntry ~= 1; err = true; end
        end
        %
        if err
          error('Wrong number of input parameters in material %s, row %d',block(1),i);
        end
        %
        params = str2double(strParams);
        switch i
          case 2
            obj.poro = params(1);
%             obj.alpha = params(2);
            obj.specGrav = params(2);
          case 3
            KTmp([1 2 4]) = [params(1), params(2), params(3)];
          case 4
            KTmp([3 5]) = [params(1), params(2)];
          case 5
            KTmp(6) = params(1);
        end
      end
      if all(KTmp([2 4 5]) == 0)
        obj.KVec = [KTmp(1); KTmp(3); KTmp(6)];
      else
        obj.KVec = KTmp;
      end
      clear KTmp
      %
      % K needs to be SPD. It is symmetric by construction but is it also
      % Positive Definite?
      K = getPermMatrix(obj);
      eigv = eig(K);
      if any(eigv < length(eigv)*eps(max(eigv)))  
        % Tolerance chosen following the hint in:
        % https://it.mathworks.com/help/matlab/math/determine-whether-matrix-is-positive-definite.html#DetermineWhetherMatrixIsSPDExample-3
        error('The permeability matrix for material %s is not positive definite',block(1));
      end
    end
  end
end