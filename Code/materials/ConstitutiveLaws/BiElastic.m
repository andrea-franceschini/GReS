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
  end

  properties (Constant, Access = private)
    STATUS_PC = 1
    STATUS_LOAD_STATE = 2

    VIRGIN = 0
    UNLOADING_RELOADING = 1
  end

  methods (Access = public)
    function obj = BiElastic(varargin)
      readMaterialParameters(obj, varargin{:});
    end

    function status = initializeStatus(obj, sigma)

      nptGauss = size(sigma, 1);
      status = zeros(nptGauss, 2);

      % Compression-negative convention: store pc as a positive pressure
      % magnitude equal to the initial mean stress magnitude.
      status(:, obj.STATUS_PC) = -mean(sigma(:, 1:3), 2);

      % The initial state lies exactly on its preconsolidation surface.
      status(:, obj.STATUS_LOAD_STATE) = obj.VIRGIN;
    end

    function [DAll, sigmaOut, statusOut] = getStiffnessMatrix( ...
        obj, sigmaIn, epsilon, dt, statusIn, cellID) %#ok<INUSD>
      % Vectorized piecewise-integrated constitutive response.
      %
      % sigmaIn and statusIn are committed quantities. epsilon is the
      % complete strain increment from the committed state to the current
      % Newton iterate.

      nptGauss = size(sigmaIn, 1);

      Kur = ratio .* obj.K;

      pcOld = statusIn(:, obj.STATUS_PC);

      % Positive scalar measures of compression.
      pOld = -mean(sigmaIn(:, 1:3), 2);
      depsVComp = -sum(epsilon(:, 1:3), 2);

      % Unloading/reloading predictor.
      pPredictor = pOld + Kur .* depsVComp;

      isVirgin = depsVComp > 0 & pPredictor > pcOld;
      
      % Default: complete increment on the unloading/reloading branch.
      pNew = pPredictor;
      Ktan = Kur;
      loadState = repmat(obj.UNLOADING_RELOADING, nptGauss, 1);

      % Points crossing into, or continuing on, virgin compression.

      % amount of reloading strain
      depsVToPc = zeros(nptGauss, 1);
      depsVToPc(isVirgin) = max( ...
        (pcOld(isVirgin) - pOld(isVirgin)) ./ Kur(isVirgin), 0);
      depsVToPc(isVirgin) = min( ...
        depsVToPc(isVirgin), depsVComp(isVirgin));

      % amount of virgin strain
      depsVVirgin = depsVComp - depsVToPc;
      pAtTransition = pOld + Kur .* depsVToPc;
      pNew(isVirgin) = pAtTransition(isVirgin) ...
        + Kvirgin(isVirgin) .* depsVVirgin(isVirgin);

      Ktan(isVirgin) = Kvirgin(isVirgin);
      loadState(isVirgin) = obj.VIRGIN;

      % Constant-deviatoric response. For engineering shear strains:
      % ds_dev,ii = 2G*(deps_ii - tr(deps)/3), ds_dev,ij = G*gamma_ij.
      depsTrace = sum(epsilon(:, 1:3), 2);
      dsDev = zeros(size(epsilon));
      dsDev(:, 1:3) = 2 .* G .* ...
        (epsilon(:, 1:3) - depsTrace ./ 3);
      dsDev(:, 4:6) = G .* epsilon(:, 4:6);

      % Increasing positive pressure means increasingly negative normal
      % stress under the compression-negative convention.
      dp = pNew - pOld;
      sigmaOut = sigmaIn + dsDev;
      sigmaOut(:, 1:3) = sigmaOut(:, 1:3) - dp;

      % Vectorized symmetric algorithmic tangent:
      % D = 2G*I_dev + Ktan*(m' m).
      DAll = zeros(6, 6, nptGauss);
      normalDiagonal = reshape(Ktan + 4 .* G ./ 3, 1, 1, []);
      normalOffDiagonal = reshape(Ktan - 2 .* G ./ 3, 1, 1, []);
      shearDiagonal = reshape(G, 1, 1, []);

      DAll(1,1,:) = normalDiagonal;
      DAll(2,2,:) = normalDiagonal;
      DAll(3,3,:) = normalDiagonal;

      DAll(1,2,:) = normalOffDiagonal;
      DAll(1,3,:) = normalOffDiagonal;
      DAll(2,1,:) = normalOffDiagonal;
      DAll(2,3,:) = normalOffDiagonal;
      DAll(3,1,:) = normalOffDiagonal;
      DAll(3,2,:) = normalOffDiagonal;

      DAll(4,4,:) = shearDiagonal;
      DAll(5,5,:) = shearDiagonal;
      DAll(6,6,:) = shearDiagonal;

      % Trial status. Commit only after global convergence.
      statusOut = statusIn;
      statusOut(:, obj.STATUS_PC) = max(pcOld, pNew);
      statusOut(:, obj.STATUS_LOAD_STATE) = loadState;

    end

    function mFactor = getMFactor(obj)
      mFactor = obj.M;
    end

    function cM = getRockCompressibility(obj, ~)

      % provisional, before refactoring the constitutive laws and status
      % management
      cM = obj.cM;

    end
  end

  methods (Access = private)
    function readMaterialParameters(obj, varargin)
      
      default = struct( ...
        'youngModulus', [], ...
        'poissonRatio', [], ...
        'bulkModulusRatio', []);

      params = readInput(default, varargin{:});

      obj.E = params.youngModulus;
      obj.nu = params.poissonRatio;
      obj.r = params.bulkModulusRatio;

      if obj.r > 1
        error('The bulk modulus ratio must be in ]0 1]')
      end

      obj.M = obj.nu ./ (1 - obj.nu);

      obj.K = obj.E ./ (3 .* (1 - 2 .* obj.nu));
      obj.G = obj.E ./ (2 .* (1 + obj.nu));
      obj.cM = 1 ./ (obj.K + 4 .* obj.G ./ 3);

    end


  end
end
