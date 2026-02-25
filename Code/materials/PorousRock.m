classdef PorousRock < handle
    % POROUS ROCK material class

    properties (Access = private)
        %General properties:
        KVec                 % Vector of permeabilities
        % (upper triangular part ordered column-wise)
        poro                 % Porosity
        biot                 % Biot coefficient
        %     alpha                % Rock compressibility (can be replaced by the oedometer test compressibility Cm)
        gamma;             % Fluid specific weight
        % specGrav             % Specific gravity of rock
        % Swr                  % Residual saturation of water
        Sr=0.;             % Residual saturation
        Ss=1.;             % Maximum saturation
    end

    methods (Access = public)
        % Class constructor method
        function obj = PorousRock(varargin)
            % Calling the function to set the object properties
            obj.readMaterialParameters(varargin{:});
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

        % function specGrav = getSpecificGravity(obj)
        %     specGrav = obj.specGrav;
        % end

        function gamma = getSpecificWeight(obj)
            gamma = obj.gamma;
        end


        % Function to get material porosity
        function biotCoeff = getBiotCoefficient(obj)
          biotCoeff = obj.biot;
        end

        % Function to get material permeability as a 3x3 matrix
        function K = getPermMatrix(obj)
          if isscalar(obj.KVec)
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
          % returns 1x9 permeability entries as
          %[Kxx,Kxy,Kxz,Kzx,Kyy,Kyz,Kxz,Kyz,Kzz]

          if isscalar(obj.KVec)
            K = [obj.KVec, 0, 0, 0, obj.KVec, 0, 0, 0, obj.KVec];
          elseif length(obj.KVec) == 3
            K = [obj.KVec(1), 0, 0, 0, obj.KVec(2), 0, 0, 0, obj.KVec(3)];
          else
            K = [obj.KVec(1), obj.KVec(6), obj.KVec(5), obj.KVec(6), obj.KVec(2), ...
              obj.KVec(4), obj.KVec(5), obj.KVec(4), obj.KVec(3)];
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
      function readMaterialParameters(obj,varargin)

        default = struct('porosity',0.3,...
                         'biotCoefficient',1.0,...
                         'permeability',[],...
                         'specificWeight',21.0,...
                         "residualSaturation",0.0,...
                         "maximumSaturation",1.0);

        % initialize also the Curve here!

        params = readInput(default,varargin{:});

        obj.poro = params.porosity;
        obj.gamma = params.specificWeight;
        obj.biot = params.biotCoefficient;
        obj.Sr = params.residualSaturation;
        obj.Ss = params.maximumSaturation;

        Kvec = params.permeability;

        nK = length(Kvec);
        if ~any([nK==1,nK==3,nK==6])
          error("Wrong number of numeric " + ...
            "values for permeability");
        end
        obj.KVec = Kvec;

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