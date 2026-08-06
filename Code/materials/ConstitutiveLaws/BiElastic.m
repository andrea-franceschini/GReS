classdef BiElastic < Elastic
  % BIELASTIC  Elastic law with two bulk moduli.
  %
  % The shear modulus is constant. The bulk modulus is:
  %   - Kvirgin during virgin compression;
  %   - Kur = r*Kvirgin during unloading/reloading.
  %
  % The transition at the preconsolidation pressure is integrated
  % piecewise. The returned stiffness is the algorithmic tangent at the
  % endpoint of the increment.
  %
  % Stress/strain convention:
  %   - Voigt ordering [xx yy zz xy yz xz];
  %   - engineering shear strains;
  %   - compression-negative stresses and strains;
  %   - p = -tr(sigma)/3 and epsVComp = -tr(epsilon) are positive in
  %     compression.
  %
  % Status layout:
  %   status(:,1) = preconsolidation pressure pc (positive magnitude)
  %   status(:,2) = loading state: 0 = virgin, 1 = unloading/reloading
  %
  % The caller passes the committed stress and status at every Newton
  % evaluation. epsilon is the current estimate of the complete strain
  % increment over the time step. The returned status is committed only
  % after global convergence.

  properties (Access = public)
    % Ratio between unloading/reloading and virgin bulk moduli:
    % Kur = r*Kvirgin. Normally r >= 1.
    r
    bulkMod
    shearMod
    DVirgin
    DUR

    nChange
  end

  properties (Constant, Access = private)
    VIRGIN = 0
    UNLOADING_RELOADING = 1
  end

  methods (Access = public)


    function obj = BiElastic(varargin)

      obj@Elastic(varargin{:});

      default = struct( ...
        'bulkModulusRatio', []);

      params = readInput(default, varargin{:});

      rat = params.bulkModulusRatio;

      if ~(rat > 0 && rat <= 1)
        error('The bulk modulus ratio must be in ]0,1]')
      end

      obj.r = 1/rat;

      obj.bulkMod = obj.E ./ (3 .* (1 - 2 .* obj.nu));
      obj.shearMod = obj.E ./ (2 .* (1 + obj.nu));

      Kur = obj.r .* obj.bulkMod;

      obj.DVirgin = obj.buildElasticMatrix(obj.bulkMod,obj.shearMod);

      obj.DUR = obj.buildElasticMatrix(Kur,obj.shearMod);


    end

    function initializeStatus(obj, varargin)

      % get stress on region local gauss points
      if nargin == 3
        domain = varargin{1};
        gaussId = varargin{2};
        stress = domain.getState.stress(gaussId,:);
      elseif nargin == 2
        stress = varargin{1};
      end

      % Compression-negative convention: store pc as a positive pressure
      % magnitude equal to the initial mean stress magnitude.
      obj.status.curr.pc = -mean(stress(:, 1:3), 2);

      % The initial state lies exactly on its preconsolidation surface.
      obj.status.curr.state = repmat(obj.VIRGIN,size(obj.status.curr.pc,1),1);

    end



    function [sigmaOut,D] = constitutiveUpdate( ...
        obj,cellList,sigmaIn,epsilon,varargin)

      gpsId = getGaussPointIds(obj.loc2gp,cellList);

      nC = numel(cellList);
      nptGauss = obj.loc2gp(cellList(1),2);
      nPoints = nptGauss.*nC;

      pcConv = obj.status.conv.pc(gpsId);

      Kvirgin = obj.bulkMod;
      G = obj.shearMod;
      Kur = obj.r.*Kvirgin;

      % Positive scalar measures of compression.
      pOld = -(1/3).*sum(sigmaIn(:,1:3),2);
      depsVComp = -sum(epsilon(:,1:3),2);

      % Unloading/reloading predictor.
      pPredictor = pOld + Kur.*depsVComp;

      tolP = 1e-8;

      isVirgin = depsVComp > 0 & ...
        pPredictor > pcConv + tolP;

      % Default: complete increment on unloading/reloading branch.
      pNew = pPredictor;

      % Amount of reloading strain required to reach pc.
      depsVToPc = zeros(nPoints,1,'like',epsilon);

      depsVToPc(isVirgin) = max( ...
        (pcConv(isVirgin)-pOld(isVirgin))./Kur,0);

      depsVToPc(isVirgin) = min( ...
        depsVToPc(isVirgin),depsVComp(isVirgin));

      % Remaining volumetric strain belongs to the virgin branch.
      depsVVirgin = depsVComp-depsVToPc;
      pAtTransition = pOld + Kur.*depsVToPc;

      pNew(isVirgin) = pAtTransition(isVirgin) + ...
        Kvirgin.*depsVVirgin(isVirgin);

      stateCurr = repmat( ...
        obj.UNLOADING_RELOADING,nPoints,1);

      stateCurr(isVirgin) = obj.VIRGIN;

      % Constant-deviatoric response.
      depsTraceOver3 = sum(epsilon(:,1:3),2)./3;
      dp = pNew-pOld;

      sigmaOut = sigmaIn;

      sigmaOut(:,1:3) = sigmaOut(:,1:3) + ...
        2.*G.*(epsilon(:,1:3)-depsTraceOver3) - dp;

      sigmaOut(:,4:6) = sigmaOut(:,4:6) + ...
        G.*epsilon(:,4:6);

      % Assign one of the two fixed tangent matrices to every material point.

      Dpages = obj.DUR(:) + ...
        (obj.DVirgin(:) - obj.DUR(:)) .* isVirgin.';

      % Preserve Gauss-point-major ordering within each cell.
      D = reshape(Dpages,6,6,nptGauss,nC);

      % Trial status. Commit only after global convergence.
      obj.status.curr.pc(gpsId) = pcConv;

      obj.status.curr.pc(gpsId(isVirgin)) = max( ...
        pcConv(isVirgin),pNew(isVirgin));

      obj.status.curr.state(gpsId) = stateCurr;

    end

    function mFactor = getMFactor(obj)
      mFactor = obj.M;
    end

    function cM = getRockCompressibility(obj, ~)

      % provisional, before refactoring the constitutive laws and status
      % management
      cM = obj.cM;

    end

    function D = buildElasticMatrix(obj,K,G)

      D = zeros(6,6);

      normalDiagonal = K + 4.*G./3;
      normalOffDiagonal = K - 2.*G./3;

      D(1:3,1:3) = normalOffDiagonal;
      D(1,1) = normalDiagonal;
      D(2,2) = normalDiagonal;
      D(3,3) = normalDiagonal;

      D(4,4) = G;
      D(5,5) = G;
      D(6,6) = G;

    end
  end


  methods (Static)

    function out = isLinear()

      out = false;

    end
  end


end
