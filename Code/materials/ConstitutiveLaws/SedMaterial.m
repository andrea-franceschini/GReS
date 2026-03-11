classdef SedMaterial < handle
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
        function obj = SedMaterial(inputStruct)
            % Calling the function to set the object properties
            obj.readMaterialParameters(inputStruct);
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
      function readMaterialParameters(obj,inputStruct)
        obj.voidRate = getXMLData(inputStruct,[],"voidRate");
        obj.gamma = getXMLData(inputStruct,[],"specificWeight");

        Kvec = getXMLData(inputStruct,[],"conductivity");
        nK = length(Kvec);
        if nK == 1
          obj.KVec(1:3) = Kvec;
        elseif nK==3
          obj.KVec = Kvec;
        else
          gresLog().error("Wrong number of numeric values for " + ...
            "hydraulic conductivity");
        end        

        obj.preStress = getXMLData(inputStruct,1.,"preStress");
        obj.inicStress = getXMLData(inputStruct,0.,"initStress");
        obj.compIdx = getXMLData(inputStruct,1,"compressibilityIndex");
        obj.rcompIdx = getXMLData(inputStruct,1,"reCompressibilityIndex");
      end
    end

    methods (Static)
      function de = getDeltaVoidRatio(sCurr,sPrev,sCons,Cc,Cr)
        % Return the variation in void ratio
          ndofs = length(sCurr);
          flag = ndofs==length(sPrev);
          flag = and(flag,ndofs==length(sCons));
          flag = and(flag,ndofs==length(Cc));
          flag = and(flag,ndofs==length(Cr));
          sCurr = abs(sCurr);
          if ~flag, return; end
          % map1 = sCurr>=sCons;
          % map2 = sCurr<=sCons;
          % map3 = and(sPrev<=sCons,map1);
          % map1 = and(map1,~map3);
          % de = zeros(ndofs,1);
          % de(map1) = -Cc(map1).*log(sCurr(map1)./sPrev(map1));
          % de(map2) = -Cr(map2).*log(sCurr(map2)./sCons(map2));
          % de(map3) = -Cc(map3).*log(sCurr(map3)./sCons(map3)) - ...
          %   Cr(map3).*log(sCons(map3)./sPrev(map3));

          map1 = sCurr <= sCons;
          map2 = sPrev >= sCons;
          map3 = and((~map1),(~map2));
          de = zeros(ndofs,1);
          de(map1) = -Cr(map1).*log10(sCurr(map1)./sPrev(map1));
          de(map2) = -Cc(map2).*log10(sCurr(map2)./sPrev(map2));          
          de(map3) = -Cr(map3).*log10(sCurr(map3)./sPrev(map3)) ...
            - Cc(map3).*log(sCurr(map3)./sCons(map3));
      end

      function de = getDevVoidRatio(sCurr,sPrev,sCons,Cc,Cr)
        % Return the variation in void ratio
          ndofs = length(sCurr);
          flag = ndofs==length(sPrev);
          flag = and(flag,ndofs==length(sCons));
          flag = and(flag,ndofs==length(Cc));
          flag = and(flag,ndofs==length(Cr));
          sCurr = abs(sCurr);
          if ~flag, return; end

          map1 = sCurr <= sCons;
          map2 = sPrev >= sCons;
          map3 = and((~map1),(~map2));
          de = zeros(ndofs,1);
          % 1/log(10) = 0.434294481903252
          de(map1) = 0.434294481903252*Cr(map1)./sCurr(map1);
          de(map2) = 0.434294481903252*Cc(map2)./sCurr(map2);
          de(map3) = 0.434294481903252*(Cr(map3)+Cc(map3))./sCurr(map3);
      end
    end
end

