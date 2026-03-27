classdef SedimentMaterial < handle
    % POROUS ROCK material class

    properties (Access = private)
        %General properties:
        compIdx            % Compressibility Index (Cc)
        rcompIdx           % Re-Compression Index (Cr)
        voidRate           % Void Rate (e0)
        preStress          % Pre Consolidate Stress(Spre)
        inicStress         % Initial stress
        KVec               % Vector of hydraulic conductivity
        gamma;             % sediment specific weight
    end

    methods (Access = public)
        % Class constructor method
        function obj = SedimentMaterial(varargin)
            % Calling the function to set the object properties
            obj.readMaterialParameters(varargin{:});
        end

        function out = getCompressibilityIdx(obj)
            out = obj.compIdx;
        end

        function out = getReCompressibilityIdx(obj)
            out = obj.rcompIdx;
        end

        function out = getVoidRate(obj)
            out = obj.voidRate;
        end

        function out = getPreConsolidadeStress(obj)
            out = obj.preStress;
        end

        function out = getInitialStress(obj)
            out = obj.inicStress;
        end

        function out = getSpecificWeight(obj)
            out = obj.gamma;
        end

        function out = getConductivity(obj)
            out = obj.KVec;
        end
    end

    methods (Access = private)
      % Assigning material parameters (check also the Materials class)
      % to object properties
      function readMaterialParameters(obj,varargin)

        % first make sure a type is defined
        default = struct('conductivity',[],...
                         'specificWeight',[],...
                         'voidRate',[],...
                         'preStress',1.,...
                         "initStress",1.,...
                         "compressibilityIndex",1.,...
                         "reCompressibilityIndex",1.);
        params = readInput(default,varargin{:});

        obj.voidRate = params.voidRate;
        obj.gamma = params.specificWeight;

        nK = length(params.conductivity);
        if nK == 1
          obj.KVec(1:3) = params.conductivity;
        elseif nK==3
          obj.KVec = params.conductivity;
        else
          gresLog().error("Wrong number of numeric values for " + ...
            "hydraulic conductivity");
        end        

        obj.preStress = params.preStress;
        obj.inicStress = params.initStress;
        obj.compIdx = params.compressibilityIndex;
        obj.rcompIdx = params.reCompressibilityIndex;
      end
    end

    methods (Static)
      % function de = getDeltaVoidRatio(sCurr,sPrev,sCons,Cc,Cr)
      %   % Return the variation in void ratio
      %     ndofs = length(sCurr);
      %     flag = ndofs==length(sPrev);
      %     flag = and(flag,ndofs==length(sCons));
      %     flag = and(flag,ndofs==length(Cc));
      %     flag = and(flag,ndofs==length(Cr));
      %     if ~flag, return; end
      %     % map = sCurr > 0; % Select only the positive stress.
      %     map = sign(sCurr) == sign(sPrev); % Select only the positive stress.
      %     map1 = and(sCurr <= sCons,map);
      %     map2 = and(sPrev >= sCons,map);
      %     map3 = and((~map1),(~map2));
      %     de = zeros(ndofs,1);
      %     de(map1) = -Cr(map1).*log10(sCurr(map1)./sPrev(map1));
      %     de(map2) = -Cc(map2).*log10(sCurr(map2)./sPrev(map2));          
      %     de(map3) = -Cr(map3).*log10(sCurr(map3)./sPrev(map3)) ...
      %       - Cc(map3).*log(sCurr(map3)./sCons(map3));
      % end

      function de = getDeltaVoidRatio(sCurr,sPrev,sCons,Cc,Cr)
        % Return the variation in void ratio
          ndofs = length(sCurr);
          flag = ndofs==length(sPrev);
          flag = and(flag,ndofs==length(sCons));
          flag = and(flag,ndofs==length(Cc));
          flag = and(flag,ndofs==length(Cr));
          if ~flag, return; end
          % map = sCurr > 0; % Select only the positive stress.
          map = sign(sCurr) == sign(sPrev); % Select only the positive stress.
          map1 = and(sCurr <= sCons,map);
          map2 = and(sPrev >= sCons,map);
          map3 = and((~map1),(~map2));
          de = zeros(ndofs,1);
          de(map1) = -Cr(map1).*log10(sCurr(map1)./sPrev(map1));
          de(map2) = -Cc(map2).*log10(sCurr(map2)./sPrev(map2));          
          de(map3) = -Cr(map3).*log10(sCons(map3)./sPrev(map3)) ...
            - Cc(map3).*log(sCurr(map3)./sCons(map3));
      end

      function de = getDevVoidRatio(sCurr,sPrev,sCons,Cc,Cr)
        % Return the variation in void ratio
          ndofs = length(sCurr);
          flag = ndofs==length(sPrev);
          flag = and(flag,ndofs==length(sCons));
          flag = and(flag,ndofs==length(Cc));
          flag = and(flag,ndofs==length(Cr));
          if ~flag, return; end
          % map = sCurr > 0; % Select only the positive stress.
          map = sign(sCurr) == sign(sPrev); % Select only the positive stress.
          map1 = and(sCurr <= sCons,map);
          map2 = and(sPrev >= sCons,map);
          map3 = and((~map1),(~map2));
          de = zeros(ndofs,1);
          % 1/log(10) = 0.434294481903252
          de(map1) = -0.434294481903252*Cr(map1)./sCurr(map1);
          de(map2) = -0.434294481903252*Cc(map2)./sCurr(map2);
          de(map3) = -0.434294481903252*Cc(map3)./sCurr(map3);
      end

      % function de = getDevVoidRatio(sCurr,sPrev,sCons,Cc,Cr)
      %   % Return the variation in void ratio
      %     ndofs = length(sCurr);
      %     flag = ndofs==length(sPrev);
      %     flag = and(flag,ndofs==length(sCons));
      %     flag = and(flag,ndofs==length(Cc));
      %     flag = and(flag,ndofs==length(Cr));
      %     if ~flag, return; end
      %     % map = sCurr > 0; % Select only the positive stress.
      %     map = sign(sCurr) == sign(sPrev); % Select only the positive stress.
      %     map1 = and(sCurr <= sCons,map);
      %     map2 = and(sPrev >= sCons,map);
      %     map3 = and((~map1),(~map2));
      %     de = zeros(ndofs,1);
      %     % 1/log(10) = 0.434294481903252
      %     de(map1) = -0.434294481903252*Cr(map1)./sCurr(map1);
      %     de(map2) = -0.434294481903252*Cc(map2)./sCurr(map2);
      %     de(map3) = -0.434294481903252*(Cr(map3)+Cc(map3))./sCurr(map3);
      % end
    end
end