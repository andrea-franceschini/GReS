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



    function [sigmaOut, D] = constitutiveUpdate( ...
        obj, cellId, sigmaIn, epsilon, varargin) 


      gpId = obj.loc2gp(cellId,1);
      nptGauss = obj.loc2gp(cellId,2);
      gpsId = gpId:gpId+nptGauss-1;
      
      pcConv = obj.status.conv.pc(gpsId);

      Kvirgin = obj.bulkMod;
      G = obj.shearMod;

      Kur = obj.r .* Kvirgin;

      % Positive scalar measures of compression.
      pOld = -(1/3)*sum(sigmaIn(:, 1:3), 2);
      depsVComp = -sum(epsilon(:, 1:3), 2);

      % Unloading/reloading predictor.
      pPredictor = pOld + Kur .* depsVComp;

      tolP = 1e-8;

      isVirgin = depsVComp > 0 & ...
        pPredictor > pcConv + tolP;



      % Default: complete increment on the unloading/reloading branch.
      pNew = pPredictor;
      Ktan = repmat(Kur,nptGauss,1);

      % Points crossing into, or continuing on, virgin compression.

      % amount of reloading strain
      depsVToPc = zeros(nptGauss, 1);

      depsVToPc(isVirgin) = max( ...
        (pcConv(isVirgin) - pOld(isVirgin)) ./ Kur, 0);
      depsVToPc(isVirgin) = min( ...
        depsVToPc(isVirgin), depsVComp(isVirgin));

      % amount of virgin strain
      depsVVirgin = depsVComp - depsVToPc;
      pAtTransition = pOld + Kur .* depsVToPc;
      pNew(isVirgin) = pAtTransition(isVirgin) ...
        + Kvirgin .* depsVVirgin(isVirgin);

      Ktan(isVirgin) = Kvirgin;
      stateCurr = repmat( ...
        obj.UNLOADING_RELOADING, numel(gpsId), 1);

      stateCurr(isVirgin) = obj.VIRGIN;

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
      D = zeros(6, 6, nptGauss);
      normalDiagonal = reshape(Ktan + 4 .* G ./ 3, 1, 1, []);
      normalOffDiagonal = reshape(Ktan - 2 .* G ./ 3, 1, 1, []);
      shearDiagonal = reshape(G, 1, 1, []);

      D(1,1,:) = normalDiagonal;
      D(2,2,:) = normalDiagonal;
      D(3,3,:) = normalDiagonal;

      D(1,2,:) = normalOffDiagonal;
      D(1,3,:) = normalOffDiagonal;
      D(2,1,:) = normalOffDiagonal;
      D(2,3,:) = normalOffDiagonal;
      D(3,1,:) = normalOffDiagonal;
      D(3,2,:) = normalOffDiagonal;

      D(4,4,:) = shearDiagonal;
      D(5,5,:) = shearDiagonal;
      D(6,6,:) = shearDiagonal;

      % Trial status. Commit only after global convergence.      
      obj.status.curr.pc(gpsId(isVirgin)) = max(pcConv, pNew);

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
  end


  methods (Static)

    function out = isLinear()

      out = false;

    end
  end


end
