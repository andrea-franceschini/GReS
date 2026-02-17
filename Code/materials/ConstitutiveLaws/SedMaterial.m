classdef SedMaterial < handle
    % POROUS ROCK material class

    properties (Access = private)
        %General properties:
        compIdx            % Compressibility index (Cc)
        rcompIdx           % Re-compression index (Cr)
        voidRate           % Void Rate (e0)
        preStress          % pre consolidate stress(Spre)
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

        

    end

    methods (Access = private)
      % Assigning material parameters (check also the Materials class)
      % to object properties
      function readMaterialParameters(obj,inputStruct)
        obj.voidRate = getXMLData(inputStruct,1,"voidRate");
        obj.preStress = getXMLData(inputStruct,1,"preStress");
        obj.compIdx = getXMLData(inputStruct,1,"compressibilityIndex");
        obj.rcompIdx = getXMLData(inputStruct,1,"reCompressibilityIndex");
      end
    end

    methods (Static)
      function out = getVoidRatio(sNew,sOld,sPre,Cc,Cr)
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
          out = zeros(ndofs,1);
          out(map1) = Cc(map1).*log(sNew(map1)./sOld(map1));
          out(map2) = Cr(map2).*log(sNew(map2)./sPre(map2));
          out(map3) = Cc(map3).*log(sNew(map3)./sPre(map3)) + ...
            Cr(map3).*log(sPre(map3)./sOld(map3));
      end
    end
end

