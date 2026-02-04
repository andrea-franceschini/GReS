classdef SedMaterial < handle
    % POROUS ROCK material class

    properties (Access = private)
        %General properties:
        compIdx            % Compressibility index (Cc)
        rcompIdx           % Re-compression index (Cr)
        voidRate           % Void Rate (e0)
        preStress          % pre consolidade stress(Spre)
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
        obj.compIdx = getXMLData(inputStruct,[],"compressibilityIndex");
        obj.rcompIdx = getXMLData(inputStruct,1,"reCompressibilityIndex");
      end
    end
end

