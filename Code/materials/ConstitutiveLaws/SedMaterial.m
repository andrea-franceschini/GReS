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
      function out = getVoidRatio(sNew,sOld,sPre,void0,Cc,Cr)
          ndofs = length(sNew);
          flag = ndofs==length(sOld);
          flag = and(flag,ndofs==length(sPre));
          flag = and(flag,ndofs==length(void0));
          flag = and(flag,ndofs==length(Cc));
          flag = and(flag,ndofs==length(Cr));
          if ~flag, return; end
          map1 = sNew>=sPre;
          map2 = ~map1;
          out = zeros(ndofs,1);
          out(map1) = void0(map1)-Cc(map1).*log(sNew(map1)./sPre(map1));
          out(map2) = void0(map2)-Cr(map2).*log(sNew(map2)./sOld(map2));
      end
    end
end

