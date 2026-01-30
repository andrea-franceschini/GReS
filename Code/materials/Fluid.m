classdef Fluid < handle
    % FLUID general fluid class

    properties (Access = private)
        %General properties:
        name
        gamma;             % Fluid specific weight
        beta;              % Fluid compressibility
        mu;                % Fluid dynamic viscosity
    end

    methods (Access = public)
      function obj = Fluid(inputStruct)
            %FLUID Class constructor method
            % Calling the function to set the material parameters
            obj.readMaterialParameters(inputStruct);
        end

        function gamma = getSpecificWeight(obj)
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
 
        function setName(obj,name)
          obj.name = name;
        end
    end

    methods (Access = private)
      function readMaterialParameters(obj,inputStruct)
            
            % Assign object properties
            obj.gamma = getXMLData(inputStruct,0,"specificWeight");
            obj.beta = getXMLData(inputStruct,0,"compressibility");
            obj.mu = getXMLData(inputStruct,[],"dynamicViscosity");
     
        end
    end
end