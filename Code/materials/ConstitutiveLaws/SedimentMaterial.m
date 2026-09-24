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

    function cb = computeOedoComp(obj,S,Sp,void,frac)
      map4 = S < obj.S1;
      map3 = (S >= obj.S1) & (S < Sp);
      map2 = (S >= Sp) & (S <= obj.S2);
      map1 = S > obj.S2;

      % map2 = S < Sp;
      % map3 = ~map2;
      
      cb = zeros(length(S),1);
      cb(map1) = obj.Cr./(log(10)*obj.S1);
      cb(map2) = obj.Cr./(log(10)*S(map2));
      cb(map3) = obj.Cc./(log(10)*S(map3));
      cb(map4) = (obj.Cc.*exp(obj.kDecay.*(obj.S2-S(map4))))./(log(10)*obj.S2);

      % voidDiff = 1+frac.*void;
      voidDiff = 1+void;
      % voidDiff(map1)=1+frac(map1).*obj.e1;
      cb=cb./voidDiff;
    end

   
    function e = getVoidRatio(obj,Scurr,Sp,ep)
      map1 = Scurr < obj.S1;
      map2 = (Scurr >= obj.S1) & (Scurr < Sp);
      map3 = (Scurr >= Sp) & (Scurr <= obj.S2);
      map4 = Scurr > obj.S2;

      e = zeros(length(Scurr),1);
      e(map1) = (1+obj.e1).*exp(obj.cbmin.*(obj.S1-Scurr(map1)))-1;
      e(map2) = ep(map2)-obj.Cr.*log10(Scurr(map2)./Sp(map2));
      e(map3) = ep(map3)-obj.Cc.*log10(Scurr(map3)./Sp(map3));
      e(map4) = obj.emin+(obj.e2-obj.emin).*exp(obj.kDecay.*(obj.S2-Scurr(map4)));
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

    function oedo = OedoCompressibility(Scurr,Sprev,Sp,void,Cc,Cr)
      % Return the variation in void ratio
      ndofs = length(Scurr);
      map1 = Scurr <= Sp;
      map2 = Sprev >= Sp;
      map3 = and((~map1),(~map2));

      de = zeros(ndofs,1);
      de(map1) = Cr(map1)./(log(10)*Scurr(map1));
      de(map2) = Cc(map2)./(log(10)*Scurr(map2));
      de(map3) = Cc(map3)./(log(10)*Scurr(map3));
      oedo = -(1./(1+void)).*de;
    end

    function oedo = OedoCompressibility2(Scurr,Sp,void,Cr,Cc)
      % Return the variation in void ratio
      ndofs = length(Scurr);
      map1 = Scurr < Sp;
      map2 = ~map1;

      de = zeros(ndofs,1);
      de(map1) = Cr(map1)./(log(10)*Scurr(map1));
      de(map2) = Cc(map2)./(log(10)*Scurr(map2));
      oedo = de./(1+void);
    end

    function dvoid = getDeltaVoidRatio(Scurr,Sprev,Sp,Cc,Cr)
      % Return the variation in void ratio
      Scurr=abs(Scurr);
      Sprev=abs(Sprev);
      Sp=abs(Sp);

      ndofs = length(Scurr);
      map1 = Scurr <= Sp;
      map2 = Sprev >= Sp;
      map3 = and((~map1),(~map2));

      dvoid = zeros(ndofs,1);
      dvoid(map1) = -Cr(map1).*log10(Scurr(map1)./Sprev(map1));
      dvoid(map2) = -Cc(map2).*log10(Scurr(map2)./Sprev(map2));
      dvoid(map3) = -Cr(map3).*log10(Sp(map3)./Sprev(map3)) ...
        - Cc(map3).*log(Scurr(map3)./Sp(map3));
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
      %     mat(1) = struct('Cc',  4,'Cr',0.4,'Sp',100,'Smin',0.1,'Smax',1e5,'emin',6,'e0',15,'S0',1);
      %     mat(2) = struct('Cc',  3,'Cr',0.1,'Sp',100,'Smin',0.1,'Smax',1e5,'emin',2.5,'e0',10,'S0',1);
      %     mat(3) = struct('Cc',0.5,'Cr',0.05,'Sp',100,'Smin',0.1,'Smax',1e5,'emin',1.5,'e0',3,'S0',1);
      %     mat(4) = struct('Cc',1.e-5,'Cr',1.e-6,'Sp',100,'Smin',0.1,'Smax',1e5,'emin',0.999975,'e0',1,'S0',1);
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



%
% mat(1) = struct('Cc',  4,'Cr',0.4,'Sp',100,'Smin',0.1,'Smax',1e4,'emin',6,'e0',15,'S0',1);
% mat(2) = struct('Cc',  3,'Cr',0.1,'Sp',100,'Smin',0.1,'Smax',1e4,'emin',2.5,'e0',10,'S0',1);
% mat(3) = struct('Cc',0.5,'Cr',0.05,'Sp',100,'Smin',0.1,'Smax',1e4,'emin',1.5,'e0',3,'S0',1);
% mat(4) = struct('Cc',1.e-5,'Cr',1.e-6,'Sp',100,'Smin',0.1,'Smax',1e4,'emin',0.999975,'e0',1,'S0',1);
% SedimentMaterial.graphVoidOedo(mat,[-3,6],1000);
%
%
%
% 
% mat(1) = struct('Cc',  4,'Cr',0.4,'Sp',100,'Smin',1e-1,'Smax',1e5,'emin',2,'e0',15,'S0',1);
% mat(2) = struct('Cc',  3,'Cr',0.1,'Sp',100,'Smin',1e-1,'Smax',1e5,'emin',1.,'e0',10,'S0',1);
% mat(3) = struct('Cc',0.5,'Cr',0.05,'Sp',100,'Smin',1e-1,'Smax',1e5,'emin',1.2,'e0',3,'S0',1);
% mat(4) = struct('Cc',1.e-5,'Cr',1.e-6,'Sp',100,'Smin',1e-1,'Smax',1e5,'emin',0.99995,'e0',1,'S0',1);
% SedimentMaterial.graphVoidOedo(mat,[-6,6],1000);