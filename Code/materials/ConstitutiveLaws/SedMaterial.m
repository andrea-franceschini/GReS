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
      function de = getVoidRatio(sNew,sOld,sPre,Cc,Cr)
        %return the variation in void ratio
          ndofs = length(sNew);
          flag = ndofs==length(sOld);
          flag = and(flag,ndofs==length(sPre));
          flag = and(flag,ndofs==length(Cc));
          flag = and(flag,ndofs==length(Cr));
          sNew = abs(sNew);
          if ~flag, return; end
          map1 = sNew>=sPre;
          map2 = sNew<=sPre;
          map3 = and(sOld<=sPre,map1);
          map1 = and(map1,~map3);
          de = zeros(ndofs,1);
          de(map1) = -Cc(map1).*log(sNew(map1)./sOld(map1));
          de(map2) = -Cr(map2).*log(sNew(map2)./sPre(map2));
          de(map3) = -Cc(map3).*log(sNew(map3)./sPre(map3)) - ...
            Cr(map3).*log(sPre(map3)./sOld(map3));
      end
    end
end

