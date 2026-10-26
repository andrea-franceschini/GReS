classdef HypoPlastic < Elastic
  % HYPOPLASTIC  Janna et al. (2012) transversely isotropic law.
  %
  % The vertical z-axis is the symmetry axis. The material parameters are
  % those introduced in Eqs. (25)-(27) of Janna et al.:
  %
  %   cM     vertical oedometric compressibility; either a positive scalar
  %          or a function handle evaluated at the signed vertical stress
  %          sigma_z;
  %   beta   Eh/Ev;
  %   nuH    Poisson ratio in the horizontal isotropy plane;
  %   nuV    vertical-to-horizontal Poisson ratio used in Eq. (21);
  %   gamma  Gh/Gv = Eh/[2*(1+nuH)*Gv];
  %   s      recompression index, i.e. the ratio between the virgin and
  %          unloading/reloading compressibilities (normally s >= 1).
  %
  % On virgin loading D = D0/cM(sigma_z). On unloading/reloading the
  % effective compressibility is cM/s and D = s*D0/cM. This class stores
  % only the dimensionless constant matrix D0.
  %
  % Stress/strain convention:
  %   - Voigt ordering [xx yy zz xy yz xz];
  %   - engineering shear strains;
  %   - compression-negative stresses and strains;
  %   - cM is evaluated with sigmaIn(:,3), without changing its sign.
  %
  % Status layout:
  %   status.curr.pc    maximum vertical compression magnitude -sigma_z;
  %   status.curr.state 0 = virgin, 1 = unloading/reloading.

  properties (Access = public)
    beta
    nuH
    nuV
    gamma
    s
    M
    D0
  end

  properties (Constant, Access = private)
    VIRGIN = 0
    UNLOADING_RELOADING = 1
  end

  properties (Access = private)
    cMInput
  end

  methods (Access = public)

    function obj = HypoPlastic(varargin)

      obj@Elastic(varargin{:});

      default = struct( ...
        'cM', [], ...
        'beta', [], ...
        'nuH', [], ...
        'nuV', [], ...
        'gamma', [], ...
        's', []);

      % readInput preserves a numeric scalar or a function handle. Thus an
      % input deck can provide, for example, cM = 1e-9 or
      % cM = @(sigmaZ) myCompressibilityLaw(sigmaZ).
      params = readInput(default,varargin{:});

      obj.cMInput = params.cM;
      obj.beta = params.beta;
      obj.nuH = params.nuH;
      obj.nuV = params.nuV;
      obj.gamma = params.gamma;
      obj.s = params.s;

      obj.checkMaterialParameters();

      obj.M = obj.nuV./(1-obj.nuH);
      obj.D0 = obj.buildDimensionlessMatrix();

    end

    function initializeStatus(obj,varargin)

      if nargin == 3
        domain = varargin{1};
        gaussId = varargin{2};
        stress = domain.getState.stress(gaussId,:);
      elseif nargin == 2
        stress = varargin{1};
      else
        error('initializeStatus requires a stress array or domain and IDs.')
      end

      % The initial state lies on the vertical preconsolidation surface.
      obj.status.curr.pc = -stress(:,3);
      obj.status.curr.state = repmat( ...
        obj.VIRGIN,size(stress,1),1);

    end

    function [sigmaOut,D] = constitutiveUpdate( ...
        obj,cellList,sigmaIn,epsilon,varargin) %#ok<INUSD>

      gpsId = getGaussPointIds(obj.loc2gp,cellList);

      nC = numel(cellList);
      nptGauss = obj.loc2gp(cellList(1),2);
      nPoints = nptGauss.*nC;

      if size(sigmaIn,1) ~= nPoints || size(epsilon,1) ~= nPoints
        error(['The stress and strain arrays must contain one row per ' ...
          'Gauss point in cellList.'])
      end

      pcConv = obj.status.conv.pc(gpsId);
      cMVirgin = obj.evaluateCompressibility(sigmaIn(:,3));

      % Base increment on the virgin branch. Because the complete matrix
      % is proportional to 1/cM, the unloading/reloading increment is
      % obtained by multiplying this vector by s.
      deltaSigmaVirgin = (epsilon*obj.D0)./cMVirgin;
      deltaSigmaUR = obj.s.*deltaSigmaVirgin;

      sigmaZCompOld = -sigmaIn(:,3);
      deltaSigmaZCompUR = -deltaSigmaUR(:,3);
      sigmaZCompPredictor = sigmaZCompOld + deltaSigmaZCompUR;

      tolSigma = 1e-8;
      isVirgin = deltaSigmaZCompUR > 0 & ...
        sigmaZCompPredictor > pcConv + tolSigma;

      % Default: the entire increment is unloading/reloading. If the
      % vertical preconsolidation stress is crossed, integrate the first
      % fraction with cM/s and the remainder with cM.
      incrementScale = repmat(obj.s,nPoints,1);
      fractionUR = zeros(nPoints,1,'like',epsilon);

      fractionUR(isVirgin) = max( ...
        (pcConv(isVirgin)-sigmaZCompOld(isVirgin))./ ...
        deltaSigmaZCompUR(isVirgin),0);

      fractionUR(isVirgin) = min(fractionUR(isVirgin),1);
      incrementScale(isVirgin) = 1 + ...
        (obj.s-1).*fractionUR(isVirgin);

      sigmaOut = sigmaIn + deltaSigmaVirgin.*incrementScale;

      % As in BiElastic, return the tangent of the branch active at the
      % endpoint: D0/cM for virgin loading and s*D0/cM otherwise.
      tangentScale = obj.s + (1-obj.s).*isVirgin;
      Dpages = obj.D0(:).*(tangentScale./cMVirgin).';
      D = reshape(Dpages,6,6,nptGauss,nC);

      % Trial status; GReS commits it only after global convergence.
      obj.status.curr.pc(gpsId) = pcConv;
      obj.status.curr.pc(gpsId(isVirgin)) = max( ...
        pcConv(isVirgin),-sigmaOut(isVirgin,3));

      stateCurr = repmat(obj.UNLOADING_RELOADING,nPoints,1);
      stateCurr(isVirgin) = obj.VIRGIN;
      obj.status.curr.state(gpsId) = stateCurr;

    end

    function mFactor = getMFactor(obj)

      mFactor = obj.M;

    end

    function cM = getRockCompressibility(obj,sigmaZ)

      if nargin < 2 || isempty(sigmaZ)
        if isnumeric(obj.cMInput)
          cM = obj.cMInput;
          return
        end
        error(['sigmaZ is required when cM is supplied as a function ' ...
          'handle.'])
      end

      cM = obj.evaluateCompressibility(sigmaZ);

    end

  end

  methods (Access = private)

    function checkMaterialParameters(obj)

      values = [obj.beta,obj.nuH,obj.nuV,obj.gamma,obj.s];

      if numel(values) ~= 5 || any(~isfinite(values)) || ...
          any(cellfun(@isempty, ...
          {obj.beta,obj.nuH,obj.nuV,obj.gamma,obj.s}))
        error('beta, nuH, nuV, gamma and s must be finite scalar values.')
      end

      if obj.beta <= 0
        error('beta = Eh/Ev must be positive.')
      end
      if obj.gamma <= 0
        error('gamma = Gh/Gv must be positive.')
      end
      if obj.s < 1
        error(['The recompression index s = cM_virgin/cM_UR must be ' ...
          'greater than or equal to one.'])
      end
      if 1-obj.nuH.^2 <= 0
        error('Thermodynamic consistency requires 1-nuH^2 > 0.')
      end
      if 1-obj.nuH-2.*obj.nuV.^2./obj.beta <= 0
        error(['Thermodynamic consistency requires ' ...
          '1-nuH-2*nuV^2/beta > 0.'])
      end

      if isnumeric(obj.cMInput)
        if ~isscalar(obj.cMInput) || ~isfinite(obj.cMInput) || ...
            obj.cMInput <= 0
          error('A numeric cM must be a finite positive scalar.')
        end
      elseif ~isa(obj.cMInput,'function_handle')
        error('cM must be a positive scalar or a function handle of sigmaZ.')
      end

    end

    function D0 = buildDimensionlessMatrix(obj)

      betaLocal = obj.beta;
      nuHLocal = obj.nuH;
      nuVLocal = obj.nuV;
      gammaLocal = obj.gamma;

      den = 1-nuHLocal.^2;
      delta = betaLocal.*(1-nuHLocal)-2.*nuVLocal.^2;

      normalDiagonal = (betaLocal-nuVLocal.^2)./den;
      normalOffDiagonal = ...
        (betaLocal.*nuHLocal+nuVLocal.^2)./den;
      normalVerticalCoupling = nuVLocal./(1-nuHLocal);
      horizontalShear = delta./(2.*den);
      verticalShear = horizontalShear./gammaLocal;

      D0 = zeros(6,6);
      D0(1,1) = normalDiagonal;
      D0(2,2) = normalDiagonal;
      D0(1,2) = normalOffDiagonal;
      D0(2,1) = normalOffDiagonal;
      D0(1,3) = normalVerticalCoupling;
      D0(3,1) = normalVerticalCoupling;
      D0(2,3) = normalVerticalCoupling;
      D0(3,2) = normalVerticalCoupling;
      D0(3,3) = 1;
      D0(4,4) = horizontalShear;
      D0(5,5) = verticalShear;
      D0(6,6) = verticalShear;

    end

    function cM = evaluateCompressibility(obj,sigmaZ)

      if isnumeric(obj.cMInput)
        cM = repmat(obj.cMInput,size(sigmaZ));
      else
        try
          cM = obj.cMInput(sigmaZ);
        catch firstError
          try
            cM = arrayfun(obj.cMInput,sigmaZ);
          catch
            rethrow(firstError)
          end
        end

        if isscalar(cM)
          cM = repmat(cM,size(sigmaZ));
        elseif numel(cM) == numel(sigmaZ)
          cM = reshape(cM,size(sigmaZ));
        else
          error('The cM function must return a scalar or one value per sigmaZ.')
        end
      end

      if ~isreal(cM) || any(~isfinite(cM(:))) || any(cM(:) <= 0)
        error('The cM law must return finite positive real values.')
      end

    end

  end

  methods (Static)

    function out = isLinear()

      % The dispatch decision cannot depend on whether this instance was
      % given a scalar or a stress-dependent cM, so use the safe choice.
      out = false;

    end

  end

end
