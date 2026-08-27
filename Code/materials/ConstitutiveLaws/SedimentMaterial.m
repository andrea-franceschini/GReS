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

        % obj.cbmin = obj.ep - obj.Cr*log10(obj.Smin/obj.Sp);
        obj.cbmin = obj.Cr./(log(10)*obj.Smin*(1+obj.e1));
        obj.kDecay = (obj.Cc)/(log(10)*obj.Smax*(obj.e2-obj.emin));
      end
    end

    methods (Static)
      % % % function de = getDeltaVoidRatio(sCurr,sPrev,sCons,Cc,Cr)
      % % %   % Return the variation in void ratio
      % % %     ndofs = length(sCurr);
      % % %     flag = ndofs==length(sPrev);
      % % %     flag = and(flag,ndofs==length(sCons));
      % % %     flag = and(flag,ndofs==length(Cc));
      % % %     flag = and(flag,ndofs==length(Cr));
      % % %     if ~flag, return; end
      % % %     sCurr = abs(sCurr);
      % % %     sPrev = abs(sPrev);
      % % %     sCons = abs(sCons);
      % % %     % map = sCurr > 0; % Select only the positive stress.
      % % %     map1 = sCurr < sCons;
      % % %     map2 = sPrev >= sCons;
      % % %     map3 = and((~map1),(~map2));
      % % %     de = zeros(ndofs,1);
      % % %     de(map1) = -Cr(map1).*log10(sCurr(map1)./sPrev(map1));
      % % %     de(map2) = -Cc(map2).*log10(sCurr(map2)./sPrev(map2));          
      % % %     de(map3) = -Cr(map3).*log10(sCons(map3)./sPrev(map3)) ...
      % % %       - Cc(map3).*log(sCurr(map3)./sCons(map3));
      % % % 
      % % %     % % % % map1 = and(sCurr <= sCons,map);
      % % %     % % % % map2 = and(sPrev >= sCons,map);
      % % %     % % % % map3 = and((~map1),(~map2));
      % % %     % % % % de = zeros(ndofs,1);
      % % %     % % % % de(map1) = -Cr(map1).*log10(sCurr(map1)./sPrev(map1));
      % % %     % % % % de(map2) = -Cc(map2).*log10(sCurr(map2)./sPrev(map2));
      % % %     % % % % de(map3) = -Cr(map3).*log10(sCons(map3)./sPrev(map3)) ...
      % % %     % % % %   - Cc(map3).*log(sCurr(map3)./sCons(map3));
      % % % end

      % % % function de = getDevVoidRatio(sCurr,sPrev,sCons,Cc,Cr)
      % % %   % Return the variation in void ratio
      % % %     ndofs = length(sCurr);
      % % %     flag = ndofs==length(sPrev);
      % % %     flag = and(flag,ndofs==length(sCons));
      % % %     flag = and(flag,ndofs==length(Cc));
      % % %     flag = and(flag,ndofs==length(Cr));
      % % %     if ~flag, return; end
      % % %     % map = sCurr > 0; % Select only the positive stress.
      % % %     map = sign(sCurr) == sign(sPrev); % Select only the positive stress.
      % % %     map1 = and(sCurr <= sCons,map);
      % % %     map2 = and(sPrev >= sCons,map);
      % % %     map3 = and((~map1),(~map2));
      % % %     de = zeros(ndofs,1);
      % % %     % 1/log(10) = 0.434294481903252
      % % %     % de(map1) = -0.434294481903252*Cr(map1)./sCurr(map1);
      % % %     % de(map2) = -0.434294481903252*Cc(map2)./sCurr(map2);
      % % %     % de(map3) = -0.434294481903252*Cc(map3)./sCurr(map3);
      % % % 
      % % % 
      % % %     de(map1) = -Cr(map1)./(log(10)*sCurr(map1));
      % % %     de(map2) = -Cc(map2)./(log(10)*sCurr(map2));
      % % %     de(map3) = -Cc(map3)./(log(10)*sCurr(map3));
      % % % end

      function void = getVoidRatioFromRef(stress,stress_Ref,void_Ref,Cc)
        % Return the variation in void ratio
        void = void_Ref - Cc.*log10(stress./stress_Ref);
      end



      function e = getVoidRatio1(Scurr,Sp,S1,S2,Cc,Cr,ep,e1,e2,emin,cbmin,kdecay)
        %GETVOIDRATIO Compute the void ratio.
        %
        %   Computes the void ratio according to the stress range:
        %               / (1+e1)*exp(cbmin*(S1-Scurr))-1,        Scurr < S1
        %              | ep-Cr*log10(Scurr/Sp),           S1 <= Scurr <= Sp
        %     e(S) =  <
        %              | ep-Cc*log10(Scurr/Sp),            Sp < Scurr <= S2
        %               \ emin+(e2-emin)*exp(k*(S2-S)),         Scurr > S2
        %
        %   where S1 and S2 define the lower and upper stress limits, Sp
        %   is the preconsolidation stress, Cc and Cr are the compression
        %   and recompression indices, and emin is the minimum void ratio.
        %
        %   INPUTS:
        %     Scurr  - Current effective stress.
        %     S1,S2  - Lower and upper stress limits.
        %     Sp     - Preconsolidation stress.
        %     Cc,Cr  - Compression and recompression indices.
        %     ep     - Previous effective void ratios.
        %     e1,e2  - Void ratios at S1 and S2.
        %     emin   - Minimum void ratio.
        %     cbmin  - Low-stress exponential coefficient.
        %     kdecay - High-stress decay coefficient.
        %
        %   OUTPUT:
        %      e     - Void ratio.

        Scurr = abs(Scurr);
        Sp    = abs(Sp);
        S1    = abs(S1);
        S2    = abs(S2);
        
        map1 = Scurr < S1;
        map4 = Scurr > S2;
        map2 = and(and(Scurr < Sp,~map1),~map4);
        map3 = and(and(Scurr > Sp,~map1),~map4);

        ndofs = length(Scurr);
        e = zeros(ndofs,1);

        e(map1) = (1+e1(map1)).*exp(cbmin(map1).*(S1(map1)-Scurr(map1)))-1;
        e(map2) = ep(map2)-Cr(map2).*log10(Scurr(map2)./Sp(map2));
        e(map3) = ep(map3)-Cc(map3).*log10(Scurr(map3)./Sp(map3));
        e(map4) = emin(map4)+(e2(map4)-emin(map4)).*exp(kdecay(map4).*(S2(map4)-Scurr(map4)));
      end

      function e = getDiffVoidRatio1(Scurr,Sp,S1,S2,Cc,Cr,e1,e2,emin,cbmin,kdecay)
        %GETVOIDRATIO Compute the void ratio.
        %
        %   Computes the derivative of the void ratio with respect to stress:
        %               / -cbmin*(1+e1)*exp(cbmin*(S1-Scurr)),   Scurr < S1
        %              | -Cr/(Scurr*log(10)),                     S1 <= Scurr <= Sp
        %     de/dS =  <
        %              | -Cc/(Scurr*log(10)),                     Sp < Scurr <= S2
        %               \ -k*(e2-emin)*exp(k*(S2-Scurr)),         Scurr > S2
        %
        %   where S1 and S2 define the lower and upper stress limits, Sp
        %   is the preconsolidation stress, Cc and Cr are the compression
        %   and recompression indices, and emin is the minimum void ratio.
        %
        %   INPUTS:
        %     Scurr  - Current effective stress.
        %     S1,S2  - Lower and upper stress limits.
        %     Sp     - Preconsolidation stress.
        %     Cc,Cr  - Compression and recompression indices.
        %     e1,e2  - Void ratios at S1 and S2.
        %     emin   - Minimum void ratio.
        %     cbmin  - Low-stress exponential coefficient.
        %     kdecay - High-stress decay coefficient.
        %
        %   OUTPUT:
        %      e     - Void ratio.

        Scurr = abs(Scurr);
        Sp    = abs(Sp);
        S1    = abs(S1);
        S2    = abs(S2);
        
        map1 = Scurr < S1;
        map4 = Scurr > S2;
        map2 = and(and(Scurr < Sp,~map1),~map4);
        map3 = and(and(Scurr > Sp,~map1),~map4);

        ndofs = length(Scurr);
        e = zeros(ndofs,1);

        e(map1) = -cbmin(map1).*(1+e1(map1)).*exp(cbmin(map1).*(S1(map1)-Scurr(map1)))-1;
        e(map2) = -Cr(map2)./(log(10)*Scurr(map2));
        e(map3) = -Cc(map3)./(log(10)*Scurr(map3));
        e(map4) = -kdecay(map4).*(e2(map4)-emin(map4)).*exp(kdecay(map4).*(S2(map4)-Scurr(map4)));
      end

      function e = getVoidRatio(Scurr,Sp,S1,S2,Cc,Cr,ep,emin)
        %GETVOIDRATIO Compute the void ratio.
        %
        %   Computes the void ratio according to the stress range:
        %               / (1+e1)*exp(cbmin*(S1-Scurr))-1,        Scurr < S1
        %              | ep-Cr*log10(Scurr/Sp),           S1 <= Scurr <= Sp
        %     e(S) =  <
        %              | ep-Cc*log10(Scurr/Sp),            Sp < Scurr <= S2
        %               \ emin+(e2-emin)*exp(k*(S2-S)),         Scurr > S2
        %
        %   where S1 and S2 define the lower and upper stress limits, Sp
        %   is the preconsolidation stress, Cc and Cr are the compression
        %   and recompression indices, and emin is the minimum void ratio.
        %
        %   INPUTS:
        %     Scurr  - Current effective stress.
        %     S1,S2  - Lower and upper stress limits.
        %     Sp     - Preconsolidation stress.
        %     Cc,Cr  - Compression and recompression indices.
        %     ep     - Previous effective void ratios.
        %     e1,e2  - Void ratios at S1 and S2.
        %     emin   - Minimum void ratio.
        %     cbmin  - Low-stress exponential coefficient.
        %     kdecay - High-stress decay coefficient.
        %
        %   OUTPUT:
        %      e     - Void ratio.

        Scurr = abs(Scurr);
        Sp    = abs(Sp);
        S1    = abs(S1);
        S2    = abs(S2);

        e1 = ep - Cr.*log10(S1./Sp);
        e2 = ep - Cc.*log10(S2./Sp);
        cbmin = Cr./(log(10)*S1.*(1+e1));
        kdecay = Cc./(log(10)*S2.*(e2-emin));
        
        map1 = Scurr < S1;
        map4 = Scurr > S2;
        map2 = and(and(Scurr < Sp,~map1),~map4);
        map3 = and(and(Scurr > Sp,~map1),~map4);

        ndofs = length(Scurr);
        e = zeros(ndofs,1);

        e(map1) = (1+e1(map1)).*exp(cbmin(map1).*(S1(map1)-Scurr(map1)))-1;
        e(map2) = ep(map2)-Cr(map2).*log10(Scurr(map2)./Sp(map2));
        e(map3) = ep(map3)-Cc(map3).*log10(Scurr(map3)./Sp(map3));
        e(map4) = emin(map4)+(e2(map4)-emin(map4)).*exp(kdecay(map4).*(S2(map4)-Scurr(map4)));
      end

      function [e,cbmin] = getDiffVoidRatio(Scurr,Sp,S1,S2,Cc,Cr,ep,emin)
        %GETVOIDRATIO Compute the void ratio.
        %
        %   Computes the derivative of the void ratio with respect to stress:
        %               / -cbmin*(1+e1)*exp(cbmin*(S1-Scurr)),   Scurr < S1
        %              | -Cr/(Scurr*log(10)),                     S1 <= Scurr <= Sp
        %     de/dS =  <
        %              | -Cc/(Scurr*log(10)),                     Sp < Scurr <= S2
        %               \ -k*(e2-emin)*exp(k*(S2-Scurr)),         Scurr > S2
        %
        %   where S1 and S2 define the lower and upper stress limits, Sp
        %   is the preconsolidation stress, Cc and Cr are the compression
        %   and recompression indices, and emin is the minimum void ratio.
        %
        %   INPUTS:
        %     Scurr  - Current effective stress.
        %     S1,S2  - Lower and upper stress limits.
        %     Sp     - Preconsolidation stress.
        %     Cc,Cr  - Compression and recompression indices.
        %     e1,e2  - Void ratios at S1 and S2.
        %     emin   - Minimum void ratio.
        %     cbmin  - Low-stress exponential coefficient.
        %     kdecay - High-stress decay coefficient.
        %
        %   OUTPUT:
        %      e     - Void ratio.

        Scurr = abs(Scurr);
        Sp    = abs(Sp);
        S1    = abs(S1);
        S2    = abs(S2);

        e1 = ep - Cr.*log10(S1./Sp);
        e2 = ep - Cc.*log10(S2./Sp);
        cbmin = Cr./(log(10)*S1.*(1+e1));
        kdecay = Cc./(log(10).*S2.*(e2-emin));
        
        map1 = Scurr < S1;
        map4 = Scurr > S2;
        map2 = and(and(Scurr < Sp,~map1),~map4);
        map3 = and(and(Scurr > Sp,~map1),~map4);

        ndofs = length(Scurr);
        e = zeros(ndofs,1);

        e(map1) = -cbmin(map1).*(1+e1(map1)).*exp(cbmin(map1).*(S1(map1)-Scurr(map1)))-1;
        e(map2) = -Cr(map2)./(log(10)*Scurr(map2));
        e(map3) = -Cc(map3)./(log(10)*Scurr(map3));
        e(map4) = -kdecay(map4).*(e2(map4)-emin(map4)).*exp(kdecay(map4).*(S2(map4)-Scurr(map4)));
      end








      function de = getDeltaVoidRatio(Scurr,Sprev,S1,S2,Sp,Cc,Cr,e1,e2,emin,cbmin,kdecay,eprev)
        %GETDELTAVOIDRATIO Compute the variation in void ratio.
        %
        %   Computes the void ratio variation according to the stress range:
        %               / (1+e1)*exp(cbmin*(S1-Scurr))-1-eprev,             Scurr < S1
        %              | -Cr*log10(Scurr/Sp),                               S1 <= Scurr <= Sp
        %     e(S) =  <  -Cr*log10(Sp/Sprev)-Cc*log10(Scurr/Sp),            Scurr> Sp and Sprev < Sp
        %              | -Cc*log10(Scurr/Sp),                               Sp < Scurr <= S2
        %               \ emin+(e2-emin)*exp(-k*(S-S2))-eprev,              Scurr > S2
        %
        %   where S1 and S2 define the lower and upper stress limits, Sp is the
        %   preconsolidation stress, Cc and Cr are the compression and
        %   recompression indices, and emin is the minimum void ratio.
        %
        %   INPUTS:
        %     Scurr  - Current effective stress.
        %     Sprev  - Previous effective stress.
        %     eprev  - Previous effective void ratios.
        %     S1,S2  - Lower and upper stress limits.
        %     Sp     - Preconsolidation stress.
        %     Cc,Cr  - Compression and recompression indices.
        %     e1,e2  - Void ratios at S1 and S2.
        %     emin   - Minimum void ratio.
        %     cbmin  - Low-stress exponential coefficient.
        %     kdecay - High-stress decay coefficient.
        %
        %   OUTPUT:
        %     de     - Variation in void ratio.

        Scurr = abs(Scurr);
        Sprev = abs(Sprev);
        Sp    = abs(Sp);
        S1    = abs(S1);
        S2    = abs(S2);

        map1 = Scurr < Sp;
        map2 = Sprev >= Sp;
        map3 = and((~map1),(~map2));
        map4 = Scurr < S1;
        map5 = Scurr > S2;

        map1 = and(map1,~map4);
        map2 = and(map2,~map4);
        map3 = and(map3,~map4);

        map1 = and(map1,~map5);
        map2 = and(map2,~map5);
        map3 = and(map3,~map5);

        ndofs = length(Scurr);
        de = zeros(ndofs,1);
        de(map1) = -Cr(map1).*log10(Scurr(map1)./Sprev(map1));
        de(map2) = -Cc(map2).*log10(Scurr(map2)./Sprev(map2));
        de(map3) = -Cr(map3).*log10(Sp(map3)./Sprev(map3)) ...
          - Cc(map3).*log(Scurr(map3)./Sp(map3));

        de(map4) = (1+e1(map4)).*exp(cbmin(map4).*(S1(map4)-Scurr(map4)))-1-eprev(map4);
        de(map5) = emin(map5)-eprev(map5)+(e2(map5)-emin(map5)).*exp(kdecay(map5).*(S2(map5)-Scurr(map5)));
      end

      function de = getDevVoidRatio(Scurr,Sprev,S1,S2,Sp,Cc,Cr,e2,emin,kdecay)
        %GETDEVVOIDRATIO2 Compute the derivative of void ratio with respect to stress.
        %
        %   Computes de/dScurr according to the stress range:
        %                / -cbmin*(1+e1)*exp(cbmin*(S1-Scurr)),             Scurr < S1
        %               | -Cr/(log(10)*Scurr),                              S1 <= Scurr < Sp
        %   de/dS =    <  -Cc*log10(Scurr/Sp),                              Scurr > Sp and Sprev < Sp
        %               | -Cc/(log(10)*Scurr),                              Sp <= Scurr <= S2
        %                \ -kdecay*(e2-emin)*exp(kdecay*(S2-Scurr)),        Scurr > S2
        %
        %   For a stress path crossing the preconsolidation stress from below,
        %   i.e. Sprev < Sp and Scurr >= Sp, the derivative corresponds to the
        %   compression branch and is therefore based on Cc.
        %
        %   INPUTS:
        %     Scurr  - Current effective stress.
        %     Sprev  - Previous effective stress.
        %     S1,S2  - Lower and upper stress limits.
        %     Sp     - Preconsolidation stress.
        %     Cc,Cr  - Compression and recompression indices.
        %     e1,e2  - Void ratios at S1 and S2.
        %     emin   - Minimum void ratio.
        %     cbmin  - Low-stress exponential coefficient.
        %     kdecay - High-stress decay coefficient.
        %
        %   OUTPUT:
        %     de     - Derivative of void ratio with respect to Scurr.

        Scurr = abs(Scurr);
        Sprev = abs(Sprev);
        Sp    = abs(Sp);
        S1    = abs(S1);
        S2    = abs(S2);

        map1 = Scurr < Sp;
        map2 = Sprev >= Sp;
        map3 = and((~map1),(~map2));
        map4 = Scurr < S1;
        map5 = Scurr > S2;

        map1 = and(map1,~map4);
        map2 = and(map2,~map4);
        map3 = and(map3,~map4);

        map1 = and(map1,~map5);
        map2 = and(map2,~map5);
        map3 = and(map3,~map5);

        ndofs = length(Scurr);
        de = zeros(ndofs,1);

        de(map1) = -Cr(map1)./(log(10)*Scurr(map1));
        de(map2) = -Cc(map2)./(log(10)*Scurr(map2));
        de(map3) = -Cc(map3)./(log(10)*Scurr(map3));

        % de(map4) = -cbmin(map4)*(1+e1(map4))*exp(cbmin(map4)*(S1(map4)-Scurr(map4)));
        % de(map4) = 0;
        de(map5) = -kdecay(map5).*(e2(map5)-emin(map5)).*exp(kdecay(map5).*(S2(map5)-Scurr(map5)));
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
        %       'Smax',1e4,'emin',2.5,'e0',10,'S0',1);
        %     mat(3) = struct('Cc',0.5,'Cr',0.05,'Sp',100,'Smin',0.1,...
        %       'Smax',1e4,'emin',1.5,'e0',3,'S0',1);
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