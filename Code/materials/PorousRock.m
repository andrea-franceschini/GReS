classdef PorousRock < handle
    % POROUS ROCK material class

    properties (Access = private)
        %General properties:
        permeability                 % Vector of permeabilities
        % (upper triangular part ordered column-wise)
        poro                 % Porosity
        biot                 % Biot coefficient
        alpha                % Rock compressibility (can be replaced by the oedometer test compressibility Cm)
        gamma;               % Fluid specific weight
        % specGrav           % Specific gravity of rock
        Sr=0.;               % Residual saturation
        Ss=1.;               % Maximum saturation        
    end


    properties
      Curves
    end

    methods (Access = public)
        % Class constructor method
        function obj = PorousRock(fluid,varargin)
            % Calling the function to set the object properties
            obj.readMaterialParameters(fluid,varargin{:});
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
          %GETSR Function to get the residual saturation of the fluid.
          Sr = obj.Sr;
        end

        function comp = getCompressibility(obj)
          comp = obj.alpha;
        end

        function gamma = getSpecificWeight(obj)
          gamma = obj.gamma;
        end


        % Function to get material porosity
        function biotCoeff = getBiotCoefficient(obj)
          biotCoeff = obj.biot;
        end

        % Function to get material permeability as a 3x3 matrix
        function K = getPermMatrix(obj,Kin)

          % convert input permeability in matrix format
          % the input should be a row vector

          if nargin == 1
            Kin = obj.permeability;
          end

          assert(size(Kin,1)==1);

          if isscalar(Kin)
            K = diag(Kin*ones(3,1));
          elseif length(Kin) == 3
            K = diag(Kin);
          elseif length(Kin) == 6
            K = [Kin(1) Kin(6) Kin(5);
              Kin(6) Kin(2) Kin(4);
              Kin(5) Kin(4) Kin(3)];
          end
        end

        function K = getPermVector(obj,Kin)
          % returns 1x9 permeability entries as
          %[Kxx,Kxy,Kxz,Kzx,Kyy,Kyz,Kxz,Kyz,Kzz]

          if nargin == 1
            Kin = obj.permeability;
          end

          nK = size(Kin,1);
          nC = size(Kin,2);
          Z = zeros(nK,1);



          if nC == 1
            K = [Kin(:), Z, Z, Z, Kin(:), Z, Z, Z, Kin(:)];
          elseif nC == 3
            K = [Kin(:,1), Z, Z, Z, Kin(:,2), Z, Z, Z, Kin(:,3)];
          elseif nC == 6
            K = [Kin(:,1), Kin(:,6), Kin(:,5), Kin(:,6), Kin(:,2), ...
              Kin(:,4), Kin(:,5), Kin(:,4), Kin(:,3)];
          end
        end

        function addCapillaryCurves(obj,f,varargin)

          fluidWeight = f.getSpecificWeight;

          obj.Curves = VanGenuchten(fluidWeight,varargin{:});
        end

        % Function to get rock compressibility
        %     function a = getRockCompressibility(obj)
        %       a = obj.alpha;
        %     end
    end

    methods (Access = private)
      % Assigning material parameters (check also the Materials class)
      % to object properties
      function readMaterialParameters(obj,fluid,varargin)

        default = struct('porosity',0.3,...
                         'biotCoefficient',1.0,...
                         'permeability',1e-12,...
                         'compressibility',0,...
                         'specificWeight',21.0,...
                         "residualSaturation",0.0,...
                         "maximumSaturation",1.0,...
                         "Curves",missing);

        % initialize also the Curve here!

        params = readInput(default,varargin{:});

        obj.poro = params.porosity;
        obj.gamma = params.specificWeight;
        obj.biot = params.biotCoefficient;
        obj.Sr = params.residualSaturation;
        obj.Ss = params.maximumSaturation;
        obj.alpha = params.compressibility;

        Kin = params.permeability;

        nK = size(Kin,2);
        if ~any([nK==1,nK==3,nK==6])
          error("Wrong number of numeric " + ...
            "values for permeability");
        end
        obj.permeability = Kin;



        % K needs to be SPD. It is symmetric by construction but is it also
        % Positive Definite?
        if size(Kin,1) == 1
          K = getPermMatrix(obj);
          eigv = eig(K);
          if any(eigv < length(eigv)*eps(max(eigv)))
            % Tolerance chosen following the hint in:
            % https://it.mathworks.com/help/matlab/math/determine-whether-matrix-is-positive-definite.html#DetermineWhetherMatrixIsSPDExample-3
            gresLog().error('The values passed for permeability matrix is not positive definite');
          end
        end

        % read capillaery curves
        if ~ismissing(params.Curves)
          addCapillaryCurves(obj,fluid,params.Curves)
        end
      end
    end
end