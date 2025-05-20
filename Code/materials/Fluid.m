classdef Fluid < handle
    % FLUID general fluid class

    properties (Access = private)
        %General properties:
        gamma;             % Fluid specific weight
        beta;              % Fluid compressibility
        mu;                % Fluid dynamic viscosity
    end

    methods (Access = public)
        function obj = Fluid(fID, matFileName)
            %FLUID Class constructor method
            % Calling the function to set the material parameters
            obj.readMaterialParameters(fID, matFileName);
        end

        function gamma = getFluidSpecWeight(obj)
            %GETFLUIDSPECWIEGHT Function to get fluid specific gravity
            gamma = obj.gamma;
        end

        function beta = getFluidCompressibility(obj)
            % GETFLUIDCOMPRESSIBILITY Function to get fluid compressibility
            beta = obj.beta;
        end

        function mu = getDynViscosity(obj)
            %GETDYNVISCOSITY Function to get fluid dynamic viscosity
            mu = obj.mu;
        end
    end

    methods (Access = private)
        function readMaterialParameters(obj, fID, matFileName)
            %READMATERIALPARAMETERS Assigning material parameters
            % (check also the Materials class) to object properties
            tmpVec = readDataInLine(fID, matFileName, 3);
            % Assign object properties
            obj.gamma = tmpVec(1);
            obj.beta = tmpVec(2);
            obj.mu = tmpVec(3);
            % if model.isVariabSatFlow()
            %     tmpVec = readDataInLine(fID, matFileName, 2);
            %     obj.Sr = tmpVec(1);
            %     obj.Ss = tmpVec(2);
            % end            
        end
    end
end