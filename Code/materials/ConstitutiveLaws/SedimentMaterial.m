classdef SedimentMaterial < handle
  % Sediment Material class

  properties (Access = public)
    %General properties:
    Cc     % Compressibility Index
    Cr     % Re-Compressibility Index
    Sp     % Pre Consolidate Stress
    S1     % Lower-stress safeguard threshold
    S2     % Upper-stress safeguard threshold
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
      obj.S1 = params.stressMin;
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
        obj.S2 = obj.Sp*10^((obj.ep - params.voidRateTrans)./obj.Cc);
      else
        obj.S2 = params.stressMax;
      end

      obj.e1 = obj.ep - obj.Cr*log10(obj.S1/obj.Sp);
      obj.e2 = obj.ep - obj.Cc*log10(obj.S2/obj.Sp);

      obj.cbmin = obj.Cr./(log(10)*obj.S1*(1+obj.e1));
      obj.kDecay = (obj.Cc)/(log(10)*obj.S2*(obj.e2-obj.emin));
    end
  end

  methods (Static)
    function out = getVoidPreCon(S,Sp,void,Cr,Cc)
      % Return the variation in void ratio
      map1 = S < Sp;
      map2 = ~map1;
      out = zeros(length(S),1);
      out(map1) = void(map1) + Cr(map1).*log10(S(map1)./Sp(map1));
      out(map2) = void(map2) + Cc(map2).*log10(S(map2)./Sp(map2));
    end

    function void = getVoidRatioFromRef(stress,stress_Ref,void_Ref,Cc)
      % Return the variation in void ratio
      void = void_Ref - Cc.*log10(stress./stress_Ref);
    end

    % % % function dvoid = getDeltaVoidRatio(Scurr,Sprev,Sp,Cc,Cr)
    % % %   % Return the variation in void ratio
    % % %   Scurr=abs(Scurr);
    % % %   Sprev=abs(Sprev);
    % % %   Sp=abs(Sp);
    % % % 
    % % %   ndofs = length(Scurr);
    % % %   map1 = Scurr < Sp;
    % % %   map2 = Sprev >= Sp;
    % % %   map3 = and((~map1),(~map2));
    % % % 
    % % %   dvoid = zeros(ndofs,1);
    % % %   dvoid(map1) = -Cr(map1).*log10(Scurr(map1)./Sprev(map1));
    % % %   dvoid(map2) = -Cc(map2).*log10(Scurr(map2)./Sprev(map2));
    % % %   dvoid(map3) = -Cr(map3).*log10(Sp(map3)./Sprev(map3)) ...
    % % %     - Cc(map3).*log(Scurr(map3)./Sp(map3));
    % % % end

    % % % function de = getDevVoidRatio(sCurr,sPrev,sCons,Cc,Cr)
    % % %   % Return the variation in void ratio
    % % %   ndofs = length(sCurr);
    % % %   flag = ndofs==length(sPrev);
    % % %   flag = and(flag,ndofs==length(sCons));
    % % %   flag = and(flag,ndofs==length(Cc));
    % % %   flag = and(flag,ndofs==length(Cr));
    % % %   if ~flag, return; end
    % % %   % map = sCurr > 0; % Select only the positive stress.
    % % %   map = sign(sCurr) == sign(sPrev); % Select only the positive stress.
    % % %   map1 = and(sCurr <= sCons,map);
    % % %   map2 = and(sPrev >= sCons,map);
    % % %   map3 = and((~map1),(~map2));
    % % % 
    % % %   de = zeros(ndofs,1);
    % % %   de(map1) = -Cr(map1)./(log(10)*sCurr(map1));
    % % %   de(map2) = -Cc(map2)./(log(10)*sCurr(map2));
    % % %   de(map3) = -Cc(map3)./(log(10)*sCurr(map3));
    % % % end

    function map = curveBranch(S,Sp,S1,S2)
      map(:,1) = S < S1;
      map(:,2) = (S >= S1) & (S <= Sp);
      map(:,3) = (S > Sp) & (S <= S2);
      map(:,4) = S > S2;
      % count = [sum(map(:,1)),sum(map(:,2)),sum(map(:,3)),sum(map(:,4))];
    end

    function void = computeVoid(S,Sp,ep,map,mat)
      void = zeros(length(S),1);
      if any(map(:,1))
        void(map(:,1)) = (1+mat.e1(map(:,1))).*exp(mat.cbmin(map(:,1)).*(mat.S1(map(:,1))-S(map(:,1))))-1;
      end
      if any(map(:,2))
        void(map(:,2)) = ep(map(:,2))-mat.Cr(map(:,2)).*log10(S(map(:,2))./Sp(map(:,2)));
      end
      if any(map(:,3))
        void(map(:,3)) = ep(map(:,3))-mat.Cc(map(:,3)).*log10(S(map(:,3))./Sp(map(:,3)));
      end
      if any(map(:,4))
        void(map(:,4)) = mat.emin(map(:,4))+(mat.e2(map(:,4))-mat.emin(map(:,4))).*exp(mat.kDecay(map(:,4)).*(mat.S2(map(:,4))-S(map(:,4))));
      end
    end

    function cb = computeOedo(S,void,map,mat)
      cb = zeros(length(S),1);
      if any(map(:,1))
        cb(map(:,1)) = mat.Cr(map(:,1))./(log(10)*mat.S1(map(:,1)));
      end
      if any(map(:,2))
        cb(map(:,2)) = mat.Cr(map(:,2))./(log(10)*S(map(:,2)));
      end
      if any(map(:,3))
        cb(map(:,3)) = mat.Cc(map(:,3))./(log(10)*S(map(:,3)));
      end
      if any(map(:,4))
        cb(map(:,4)) = (mat.Cc(map(:,4)).*exp(mat.kDecay(map(:,4)).*(mat.S2(map(:,4))-S(map(:,4)))))./(log(10)*mat.S2(map(:,4)));
      end

      voidDiff = 1+void;
      if any(map(:,1))
        voidDiff(map(:,1))=1+mat.e1(map(:,1));
      end
      cb=cb./voidDiff;
    end

    function graphVoidOedo(mat,range,pts)
      % GRAPHVOIDOEDO Plot void ratio and oedometric compressibility.
      %
      %   graphVoidOedo(mat,range,pts) computes the void ratio and oedometric
      %   compressibility for one or more sediment materials and plots the
      %   corresponding curves.
      %
      %   INPUTS:
      %     mat   - Material parameters:
      %             name optional material name used in the legend
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
      %   If the field 'name' is not provided or is empty, the default
      %   material names are 'mat1', 'mat2', ..., 'matN'.
      %
      %   EXAMPLE:
      %     mat(1) = struct('name','peat', ...
      %                     'Cc',4,'Cr',0.4,'Sp',100, ...
      %                     'Smin',0.1,'Smax',1e5, ...
      %                     'emin',6,'e0',15,'S0',1);
      %
      %     mat(2) = struct('name','clay', ...
      %                     'Cc',3,'Cr',0.1,'Sp',100, ...
      %                     'Smin',0.1,'Smax',1e5, ...
      %                     'emin',2.5,'e0',10,'S0',1);
      %
      %     mat(3) = struct('name','silt', ...
      %                     'Cc',0.5,'Cr',0.05,'Sp',100, ...
      %                     'Smin',0.1,'Smax',1e5, ...
      %                     'emin',1.5,'e0',3,'S0',1);
      %
      %     SedimentMaterial.graphVoidOedo(mat,[-3,6],1000);
      nmat = length(mat);

      % Stress values
      sigma = logspace(range(1),range(2),pts);

      % Preallocate arrays
      oedo = zeros(pts,nmat);
      void = zeros(pts,nmat);

      % Legend names
      labels = cell(1,nmat);

      for loop = 1:nmat

        % Use material name when available, otherwise use default name
        if isfield(mat,'name') && ~isempty(mat(loop).name)
          labels{loop} = mat(loop).name;
        else
          labels{loop} = sprintf('mat%d',loop);
        end

        % Void ratio at preconsolidation stress
        if mat(loop).S0 < mat(loop).Sp
          mat(loop).ep = mat(loop).e0 ...
            - mat(loop).Cr*log10(mat(loop).S0/mat(loop).Sp);
        else
          mat(loop).ep = mat(loop).e0 ...
            - mat(loop).Cc*log10(mat(loop).S0/mat(loop).Sp);
        end

        % Void ratio at stress limits
        mat(loop).e1 = mat(loop).ep ...
          - mat(loop).Cr*log10(mat(loop).Smin/mat(loop).Sp);

        mat(loop).e2 = mat(loop).ep ...
          - mat(loop).Cc*log10(mat(loop).Smax/mat(loop).Sp);

        % Low-stress exponential coefficient
        mat(loop).cbmin = mat(loop).Cr / ...
          (log(10)*mat(loop).Smin*(1+mat(loop).e1));

        % High-stress exponential decay coefficient
        mat(loop).kDecay = mat(loop).Cc / ...
          (log(10)*mat(loop).Smax* ...
          (mat(loop).e2-mat(loop).emin));
      end

      % Compute void ratio and oedometric compressibility
      for loop = 1:nmat
        for pt = 1:pts
          if sigma(pt) < mat(loop).Sp
            % Low-stress region
            if sigma(pt) < mat(loop).Smin
              void(pt,loop) = (1+mat(loop).e1) * ...
                exp(mat(loop).cbmin * (mat(loop).Smin-sigma(pt))) - 1;

              oedo(pt,loop) = mat(loop).cbmin;
            else
              void(pt,loop) = mat(loop).ep ...
                - mat(loop).Cr*log10(sigma(pt)/mat(loop).Sp);

              oedo(pt,loop) = ...
                mat(loop).Cr / (log(10)*sigma(pt)*(1+void(pt,loop)));
            end
          else
            % High-stress exponential region
            if sigma(pt) > mat(loop).Smax
              void(pt,loop) = ...
                mat(loop).emin + (mat(loop).e2-mat(loop).emin) * ...
                exp(mat(loop).kDecay*(mat(loop).Smax-sigma(pt)));

              oedo(pt,loop) = mat(loop).kDecay * ...
                (void(pt,loop)-mat(loop).emin) / (1+void(pt,loop));
            else
              void(pt,loop) = mat(loop).ep ...
                - mat(loop).Cc*log10(sigma(pt)/mat(loop).Sp);

              oedo(pt,loop) = ...
                mat(loop).Cc / (log(10)*sigma(pt)*(1+void(pt,loop)));
            end
          end
        end
      end

      % Oedometric compressibility
      figure('Position',[100,100,700,700]);
      hold on;
      h = plot(sigma,oedo,'-','LineWidth',2,'MarkerSize',14);
      for loop = 1:nmat
        xline(mat(loop).Smin,'--','\sigma_{min}','HandleVisibility','off');
        xline(mat(loop).Smax,'--','\sigma_{max}','HandleVisibility','off');
      end
      xlabel('Stress');
      ylabel('OedoComp');
      legend(h,labels,'Location','best');

      set(gca, 'FontName','Liberation Serif', ...
        'FontSize',16,'XGrid','on','YGrid','on','XScale','log');

      % Void ratio
      figure('Position',[100,100,700,700]);
      hold on;
      h = plot(sigma,void,'-','LineWidth',2,'MarkerSize',14);
      for loop = 1:nmat
        xline(mat(loop).Smin,'--','\sigma_{min}','HandleVisibility','off');
        xline(mat(loop).Smax,'--','\sigma_{max}','HandleVisibility','off');
      end
      xlabel('Stress');
      ylabel('Void Ratio');

      legend(h,labels,'Location','best');
      set(gca, 'FontName','Liberation Serif', ...
        'FontSize',16, 'XGrid','on', 'YGrid','on', 'XScale','log');
    end

  end
end



%
% mat(1) = struct('name','peat','Cc',  4,'Cr',0.4,'Sp',100,'Smin',0.1,'Smax',1e4,'emin',6,'e0',15,'S0',1);
% mat(2) = struct('name','clay','Cc',  3,'Cr',0.1,'Sp',100,'Smin',0.1,'Smax',1e4,'emin',2.5,'e0',10,'S0',1);
% mat(3) = struct('name','silt','Cc',0.5,'Cr',0.05,'Sp',100,'Smin',0.1,'Smax',1e4,'emin',1.5,'e0',3,'S0',1);
% mat(4) = struct('name','base','Cc',1.e-5,'Cr',1.e-6,'Sp',100,'Smin',0.1,'Smax',1e4,'emin',0.999975,'e0',1,'S0',1);
% SedimentMaterial.graphVoidOedo(mat,[-3,6],1000);
%
%
%
%
% mat(1) = struct('name','peat','Cc',  4,'Cr',0.4,'Sp',100,'Smin',1e-1,'Smax',1e5,'emin',2,'e0',15,'S0',1);
% mat(2) = struct('name','clay','Cc',  3,'Cr',0.1,'Sp',100,'Smin',1e-1,'Smax',1e5,'emin',1.,'e0',10,'S0',1);
% mat(3) = struct('name','silt','Cc',0.5,'Cr',0.05,'Sp',100,'Smin',1e-1,'Smax',1e5,'emin',1.2,'e0',3,'S0',1);
% mat(4) = struct('name','base','Cc',1.e-5,'Cr',1.e-6,'Sp',100,'Smin',1e-1,'Smax',1e5,'emin',0.99995,'e0',1,'S0',1);
% SedimentMaterial.graphVoidOedo(mat,[-2,6],1000);

% mat(1) = struct('name','peat','Cc',  4,'Cr',0.4,'Sp',100,'Smin',1e-1,'Smax',1e5,'emin',2,'e0',15,'S0',1);
% mat(2) = struct('name','clay','Cc',  3,'Cr',0.1,'Sp',100,'Smin',1e-1,'Smax',1e5,'emin',1.,'e0',10,'S0',1);
% mat(3) = struct('name','silt','Cc',0.5,'Cr',0.05,'Sp',100,'Smin',1e-1,'Smax',1e5,'emin',1.2,'e0',3,'S0',1);
% mat(4) = struct('name','base','Cc',1.e-5,'Cr',1.e-6,'Sp',100,'Smin',1e-1,'Smax',1e5,'emin',0.99995,'e0',1,'S0',1);
% SedimentMaterial.graphVoidOedo(mat,[-2,6],1000);