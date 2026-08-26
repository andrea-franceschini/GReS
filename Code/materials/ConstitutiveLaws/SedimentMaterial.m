classdef SedimentMaterial < handle
    % Sediment Material class

    properties (Access = public)
        %General properties:
        Cc     % Compressibility Index
        Cr     % Re-Compressibility Index
        Sp     % Pre Consolidate Stress
        Smin   % Lower-stress safeguard threshold
        Smax   % Upper-stress safeguard threshold
        emin   % Asymptotic minimum void ratio

        e0     % Reference void ratio state
        S0     % Reference stress state

        gamma  % Specific weight
        KVec   % Hydraulic conductivity
    end

    properties (Access = public)
        %General properties:
        ep     % pre-consolidated void ratio
        e1     % void ratio reference to Smin
        e2     % void ratio reference to Smax
        cbmin  % oedometric compressibility min
        kDecay % decay rate coefficient
    end

    methods (Access = public)
        % Class constructor method
        function obj = SedimentMaterial(varargin)
            % Calling the function to set the object properties
            obj.readMaterialParameters(varargin{:});
        end

        function out = getCompressibilityIdx(obj)
            out = obj.Cc;
        end

        function out = getReCompressibilityIdx(obj)
            out = obj.Cr;
        end

        function out = getSpecificWeight(obj)
            out = obj.gamma;
        end

        function out = getConductivity(obj)
            out = obj.KVec;
        end



        % % function out = getVoidRate(obj)
        % %     out = obj.voidRate;
        % % end
        % % 
        % % function out = getVoidLowerLim(obj)
        % %     out = obj.voidLowerLim;
        % % end
        % % 
        % % function out = getPreConsolidadeStress(obj)
        % %     out = obj.preStress;
        % % end
        % % 
        % % function out = getInitialStress(obj)
        % %     out = obj.inicStress;
        % % end

        
    end

    methods (Access = private)
      % Assigning material parameters (check also the Materials class)
      % to object properties
      function readMaterialParameters(obj,varargin)
        % first make sure a type is defined
        default = struct('conductivity',[],...
                         'specificWeight',[],...
                         'stressReference',[],...
                         "voidRateReference",[],...
                         "compressibilityIndex",[],...
                         "reCompressibilityIndex",[],...
                         'stressPreConsolidated',[],...
                         'stressMin',[],...
                         'stressMax',[],...
                         'voidRateMin',[]);
        params = readInput(default,varargin{:});
        if any(structfun(@isempty, params))
          gresLog().error("Material not well defined, at least" + ...
            " one field missed.");
        end

        obj.gamma = params.specificWeight;
        nK = length(params.conductivity);
        if nK == 1
          obj.KVec(1:3) = params.conductivity;
        elseif nK==3
          obj.KVec = params.conductivity;
        else
          gresLog().error("Wrong number of numeric values for " + ...
            "hydraulic conductivity");
        end

        obj.Cc = params.compressibilityIndex;
        obj.Cr = params.reCompressibilityIndex;
        obj.Sp = params.stressPreConsolidated;
        obj.Smin = params.stressMin;        
        obj.emin = params.voidRateMin;

        obj.e0 = params.voidRateReference;
        obj.S0 = params.stressReference;

        if obj.S0<obj.Sp
          obj.ep = obj.e0 - obj.Cr*log10(obj.S0/obj.Sp);
        else
          obj.ep = obj.e0 - obj.Cc*log10(obj.S0/obj.Sp);
        end

        if isfield(params,"voidRateTrans")
          if obj.ep<params.voidRateTrans
            gresLog().error("Material not well defined, void ratio for" + ...
              " the pre-consolidated structure is lower than the" + ...
              " voidRateTrans value.");
          end
          obj.Smax = obj.Sp*10^((obj.ep - params.voidRateTrans)./obj.Cc);
        else
          obj.Smax = params.stressMax;
        end

        obj.e1 = obj.ep - obj.Cr*log10(obj.Smin/obj.Sp);
        obj.e2 = obj.ep - obj.Cc*log10(obj.Smax/obj.Sp);

        obj.cbmin = obj.ep - obj.Cr*log10(obj.Smin/obj.Sp);
        obj.kDecay = (obj.Cc)/(log(10)*obj.Smax*(obj.e2-obj.emin));
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
          if ~flag, return; end
          sCurr = abs(sCurr);
          sPrev = abs(sPrev);
          sCons = abs(sCons);
          % map = sCurr > 0; % Select only the positive stress.
          map1 = sCurr < sCons;
          map2 = sPrev >= sCons;
          map3 = and((~map1),(~map2));
          de = zeros(ndofs,1);
          de(map1) = -Cr(map1).*log10(sCurr(map1)./sPrev(map1));
          de(map2) = -Cc(map2).*log10(sCurr(map2)./sPrev(map2));          
          de(map3) = -Cr(map3).*log10(sCons(map3)./sPrev(map3)) ...
            - Cc(map3).*log(sCurr(map3)./sCons(map3));

          % % % % map1 = and(sCurr <= sCons,map);
          % % % % map2 = and(sPrev >= sCons,map);
          % % % % map3 = and((~map1),(~map2));
          % % % % de = zeros(ndofs,1);
          % % % % de(map1) = -Cr(map1).*log10(sCurr(map1)./sPrev(map1));
          % % % % de(map2) = -Cc(map2).*log10(sCurr(map2)./sPrev(map2));
          % % % % de(map3) = -Cr(map3).*log10(sCons(map3)./sPrev(map3)) ...
          % % % %   - Cc(map3).*log(sCurr(map3)./sCons(map3));
      end

      function de = getDevVoidRatio(sCurr,sPrev,sCons,Cc,Cr)
        % Return the variation in void ratio
          ndofs = length(sCurr);
          flag = ndofs==length(sPrev);
          flag = and(flag,ndofs==length(sCons));
          flag = and(flag,ndofs==length(Cc));
          flag = and(flag,ndofs==length(Cr));
          if ~flag, return; end
          % map = sCurr > 0; % Select only the positive stress.
          map = sign(sCurr) == sign(sPrev); % Select only the positive stress.
          map1 = and(sCurr <= sCons,map);
          map2 = and(sPrev >= sCons,map);
          map3 = and((~map1),(~map2));
          de = zeros(ndofs,1);
          % 1/log(10) = 0.434294481903252
          % de(map1) = -0.434294481903252*Cr(map1)./sCurr(map1);
          % de(map2) = -0.434294481903252*Cc(map2)./sCurr(map2);
          % de(map3) = -0.434294481903252*Cc(map3)./sCurr(map3);


          de(map1) = -Cr(map1)./(log(10)*sCurr(map1));
          de(map2) = -Cc(map2)./(log(10)*sCurr(map2));
          de(map3) = -Cc(map3)./(log(10)*sCurr(map3));
      end

      function void = getVoidRatioFromRef(stress,stress_Ref,void_Ref,Cc)
        % Return the variation in void ratio
        void = void_Ref - Cc.*log10(stress./stress_Ref);
      end




      function graphVoidOedo(mat,range,pts)
        % GRAPHVOIDOEDO Plot void ratio and oedometric compressibility.
        % 
        %   graphVoidOedo(mat,range,pts) computes the void ratio and oedometric
        %   compressibility for the material
        % 
        %   INPUTS:
        %     mat   - Material parameters:
        %             Cc   compression index
        %             Cr   recompression index
        %             Sp   preconsolidation stress
        %             Smin minimum stress
        %             Smax maximum stress
        %             emin minimum void ratio
        %             e0   reference void ratio
        %             S0   reference stress
        % 
        %     range - Exponents defining the stress range used by LOGSPACE.
        %     pts   - Number of stress points.
        % 
        %   EXAMPLE:
        %     mat(1) = struct('Cc',  4,'Cr',0.4,'Sp',100,'Smin',0.1,...
        %       'Smax',1e4,'emin',6,'e0',15,'S0',1);
        %     mat(2) = struct('Cc',  3,'Cr',0.1,'Sp',100,'Smin',0.1,...
        %       'Smax',1e4,'emin',3,'e0',10,'S0',1);
        %     mat(3) = struct('Cc',0.5,'Cr',0.05,'Sp',100,'Smin',0.1,...
        %       'Smax',1e4,'emin',1.8,'e0',3,'S0',1);
        %     mat(4) = struct('Cc',1.e-5,'Cr',1.e-6,'Sp',100,'Smin',0.1,...
        %       'Smax',1e4,'emin',0.999975,'e0',1,'S0',1);
        %     SedimentMaterial.graphVoidOedo(mat,[-3,6],1000);

        nmat = length(mat);
        sigma = logspace(range(1),range(2),pts);
        oedo = zeros(pts,nmat);
        void = zeros(pts,nmat);
        for loop=1:nmat
          if mat(loop).S0<mat(loop).Sp
            mat(loop).ep = mat(loop).e0 - mat(loop).Cr*log10(mat(loop).S0/mat(loop).Sp);
          else
            mat(loop).ep = mat(loop).e0 - mat(loop).Cc*log10(mat(loop).S0/mat(loop).Sp);
          end
          mat(loop).e1 = mat(loop).ep - mat(loop).Cr*log10(mat(loop).Smin/mat(loop).Sp);
          mat(loop).e2 = mat(loop).ep - mat(loop).Cc*log10(mat(loop).Smax/mat(loop).Sp);

          mat(loop).cbmin = mat(loop).Cr/(log(10)*mat(loop).Smin*(1+mat(loop).e1));
          mat(loop).kDecay = mat(loop).Cc/(log(10)*mat(loop).Smax*(mat(loop).e2-mat(loop).emin));
        end

        for loop=1:nmat
          for pt = 1:pts
            if sigma(pt) < mat(loop).Sp
              if sigma(pt) < mat(loop).Smin
                void(pt,loop) = (1+mat(loop).e1)*exp(mat(loop).cbmin*(mat(loop).Smin-sigma(pt)))-1;
                oedo(pt,loop) = mat(loop).cbmin;
              else
                void(pt,loop) = mat(loop).ep - mat(loop).Cr*log10(sigma(pt)/mat(loop).Sp);
                oedo(pt,loop) = mat(loop).Cr/(log(10)*sigma(pt)*(1+void(pt,loop)));
              end
            else
              if sigma(pt) > mat(loop).Smax
                void(pt,loop) = mat(loop).emin+(mat(loop).e2-mat(loop).emin)*exp(mat(loop).kDecay*(mat(loop).Smax-sigma(pt)));
                oedo(pt,loop) = mat(loop).kDecay*(void(pt,loop)-mat(loop).emin)/(1+void(pt,loop));
              else
                void(pt,loop) = mat(loop).ep - mat(loop).Cc*log10(sigma(pt)/mat(loop).Sp);
                oedo(pt,loop) = mat(loop).Cc/(log(10)*sigma(pt)*(1+void(pt,loop)));
              end
            end
          end
        end

        figure('Position', [100, 100, 700, 700]);
        hold on;
        plot(sigma,oedo,'-', 'LineWidth', 2, 'MarkerSize', 14);
        for loop=1:nmat
          xline(mat(loop).Smin , '--', '\sigma_{min}');
          xline(mat(loop).Smax , '--', '\sigma_{max}');
        end
        xlabel('Stress');
        ylabel('OedoComp');
        set(gca,'FontName','Liberation Serif','FontSize',16,...
          'XGrid','on','YGrid','on','XScale','log');

        figure('Position', [100, 100, 700, 700]);
        hold on;
        plot(sigma,void,'-', 'LineWidth', 2, 'MarkerSize', 14);
        for loop=1:nmat
          xline(mat(loop).Smin , '--', '\sigma_{min}');
          xline(mat(loop).Smax , '--', '\sigma_{max}');
        end
        xlabel('Stress');
        ylabel('Void Rate');
        set(gca,'FontName','Liberation Serif','FontSize',16,...
          'XGrid','on','YGrid','on','XScale','log');
      end

    end
end