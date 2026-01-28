classdef PorousRock < handle
    % POROUS ROCK material class

    properties (Access = private)
        %General properties:
        KVec                 % Vector of permeabilities
        % (upper triangular part ordered column-wise)
        poro                 % Porosity
        biot                 % Biot coefficient
        %     alpha                % Rock compressibility (can be replaced by the oedometer test compressibility Cm)
        specGrav             % Specific gravity of rock
        % Swr                  % Residual saturation of water
        Sr=0.;             % Residual saturation
        Ss=1.;             % Maximum saturation
    end

    methods (Access = public)
        % Class constructor method
        function obj = PorousRock(inputStruct)
            % Calling the function to set the object properties
            obj.readMaterialParameters(inputStruct);
        end

        % Function to get material porosity
        function poro = getPorosity(obj)
            poro = obj.poro;
        end

        % Function to get material porosity
        % function Swr = getWaterResSat(obj)
        %     Swr = obj.Swr;
        % end

        function Ss = getMaxSaturation(obj)
            %GETSS Function to get the maximun saturation of the fluid.
            Ss = obj.Ss;
        end

        function Sr = getResidualSaturation(obj)
            %GETSS Function to get the residual saturation of the fluid.
            Sr = obj.Sr;
        end

        function specGrav = getSpecificGravity(obj)
            specGrav = obj.specGrav;
        end


        % Function to get material porosity
        function biotCoeff = getBiotCoefficient(obj)
            biotCoeff = obj.biot;
        end

        % Function to get material permeability as a 3x3 matrix
        function K = getPermMatrix(obj)
            if length(obj.KVec) == 1
                K = diag(obj.KVec*ones(3,1));
            elseif length(obj.KVec) == 3
                K = diag(obj.KVec);
            else
                K = [obj.KVec(1) obj.KVec(6) obj.KVec(5);
                    obj.KVec(6) obj.KVec(2) obj.KVec(4);
                    obj.KVec(5) obj.KVec(4) obj.KVec(3)];
            end
        end

        function K = getPermVector(obj)
            if isscalar(obj.KVec)  % output [K K K 0 0 0]
              K = obj.KVec*[ones(1,3) zeros(1,3)];
            elseif length(obj.KVec) == 3  % output [Kxx Kyy Kzz 0 0 0]
              K = [obj.KVec zeros(1,3)];
            else  % output [Kxx Kyy Kzz Kyz Kxz Kxy]
              K = obj.KVec;
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
      function readMaterialParameters(obj,inputStruct)

        obj.poro = getXMLData(inputStruct,[],"porosity");
        obj.specGrav = getXMLData(inputStruct,21,"specificGravity");
        obj.biot = getXMLData(inputStruct,1,"biotCoefficient");
        Kvec = getXMLData(inputStruct,[],"permeability");
        nK = length(Kvec);
        if ~any([nK==1,nK==3,nK==6])
          gres_log().error("Wrong number of numeric " + ...
            "values for permeability");
        end
        obj.KVec = Kvec;

        obj.Sr = getXMLData(inputStruct,0,"residualSaturation");
        obj.Ss = getXMLData(inputStruct,1,"maximumSaturation");

        % K needs to be SPD. It is symmetric by construction but is it also
        % Positive Definite?
        K = getPermMatrix(obj);
        eigv = eig(K);
        if any(eigv < length(eigv)*eps(max(eigv)))
          % Tolerance chosen following the hint in:
          % https://it.mathworks.com/help/matlab/math/determine-whether-matrix-is-positive-definite.html#DetermineWhetherMatrixIsSPDExample-3
          error('The permeability matrix for material %s is not positive definite',matFileName);
        end
      end
    end
end