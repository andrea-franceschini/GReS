classdef DruckerPrager < handle
  % DRUCKERPRAGER GEOS-style Drucker-Prager elastoplastic material.
  %
  % The yield function is written as
  %
  %   F = q + alpha*p - cohesion <= 0,
  %
  % where q = sqrt(3*J2), p = tr(sigma)/3, alpha is the friction
  % coefficient, beta is the dilation coefficient, and cohesion is the
  % current cone intercept. 
  % The internal hardening variable stored in the first column
  % of the status array is the current cohesion.

  properties (Access = public)
    % Elastic modulus
    E
    % Poisson's ratio
    nu
    % M factor (sigmax/sigmaz)
    M
    % Vertical compressibility Cm
    cM
    % Dilatancy angle
    psi
    % Friction angle
    phi
    % Input Mohr-Coulomb cohesion c
    co
    % Hardening parameter acting directly on the GEOS cohesion/intercept
    h
    % Drucker-Prager coefficients
    alpha       % friction coefficient in F = q + alpha*p - cohesion
    beta        % dilation coefficient in G = q + beta*p
    epsilon     % conversion from Mohr-Coulomb cohesion c to DP intercept
  end

  methods (Access = public)
    function obj = DruckerPrager(varargin)
      obj.readMaterialParameters(varargin{:});
    end

    function cohesion = initializeStatus(obj, sigma)
      nptGauss = size(sigma,1);
      cohesion = zeros(nptGauss,6);
      cohesion(:,1) = obj.epsilon * obj.co;
    end

    function [DAll, sigmaOut, cohesion] = getStiffnessMatrix(obj, sigmaIn, epsilonIn, dt, convCohesion, cellID) %#ok<INUSD>
      
      nptGauss = size(sigmaIn,1);

      % if isempty(cohesion) || isscalar(cohesion)
      %   cohesion = obj.initializeStatus(sigmaIn);
      % end

      D = obj.getElasticTensor();
      DAll = repmat(D,[1, 1, nptGauss]);

      % Incoming epsilonIn is an engineering strain increment row matrix.
      sigmaTrial = sigmaIn + epsilonIn * D;
      [sigmaOut, DAll, cohesion] = obj.returnMapping(sigmaTrial, nptGauss, DAll, convCohesion);
    end

    function m = getMFactor(obj)
      m = obj.M;
    end

    function cM = getRockCompressibility(obj)
      cM = obj.cM;
    end




    function [sigma, D, oldCohesion] = returnMapping(obj, sigmaTrial, nptGauss, D, oldCohesion)

      sigma = sigmaTrial;

      I = [1;1;1;0;0;0];
      Pdev = diag([1, 1, 1, 0.5, 0.5, 0.5]) - (1/3) * (I * I');

      G = obj.E / (2 * (1 + obj.nu));
      K = obj.E / (3 * (1 - 2 * obj.nu));

      tolYield = 1e-9;
      tolQ     = 1e-12;
      tolAlpha = 1e-14;

      for i = 1:nptGauss

        pTrial = sum(sigmaTrial(i,1:3)) / 3;
        sTrial = sigmaTrial(i,1:6)' - pTrial * I;

        qTrial = sqrt(0.5 * ((sigmaTrial(i,1) - sigmaTrial(i,2))^2 + ...
          (sigmaTrial(i,1) - sigmaTrial(i,3))^2 + ...
          (sigmaTrial(i,2) - sigmaTrial(i,3))^2) + ...
          3.0 * (sigmaTrial(i,4)^2 + sigmaTrial(i,5)^2 + sigmaTrial(i,6)^2));

        cOld = oldCohesion(i,1);

        yield = qTrial + obj.alpha * pTrial - cOld;

        if yield < tolYield
          continue;
        end

        % Standard GEOS-style cone return.
        denom = 3.0 * G + obj.alpha * obj.beta * K + obj.h;
        lambdaCone = yield / denom;

        cCone = cOld + lambdaCone * obj.h;
        if cCone < 0
          cCone = 0;
        end

        pCone = pTrial - lambdaCone * K * obj.beta;
        qCone = qTrial - lambdaCone * 3.0 * G;

        % The regular cone return is valid only if q remains positive.
        useApex = (qTrial <= tolQ) || (qCone <= tolQ);

        if ~useApex

          % smooth return

          sigma(i,1:6) = (pCone * I + (qCone / qTrial) * sTrial)';

          theta = obj.alpha * K * I + (3.0 * G / qTrial) * sTrial;
          ppsi  = obj.beta  * K * I + (3.0 * G / qTrial) * sTrial;

          D(:,:,i) = K * (I * I') + ...
            2.0 * G * (qCone / qTrial) * Pdev + ...
            (3.0 * G / qTrial^2) * (1.0 - qCone / qTrial) * (sTrial * sTrial') - ...
            (1.0 / denom) * (ppsi * theta');

          oldCohesion(i,1) = cCone;

        else


          % Apex return

          if abs(obj.alpha) < tolAlpha
            % Degenerate case: alpha = 0 gives a von-Mises-like cylinder,
            % so there is no Drucker-Prager apex. Fall back to q = 0 only.
            pApex = pTrial;
            cNew  = cOld;
            D(:,:,i) = K * (I * I');

          else

            denomApex = obj.alpha * obj.beta * K + obj.h;

            if abs(denomApex) > tolAlpha
              lambdaApex = (obj.alpha * pTrial - cOld) / denomApex;

              if lambdaApex < 0
                lambdaApex = 0;
              end

              cNew = cOld + lambdaApex * obj.h;

              if cNew < 0
                cNew = 0;
              end

              pApex = cNew / obj.alpha;

              % Consistent tangent for the active apex branch.
              %
              D(:,:,i) = (K * obj.h / denomApex) * (I * I');

            else
              % Perfectly plastic, no volumetric plastic flow / no hardening.
              % The apex equation is singular. Use a robust projection to the
              % current apex and a zero tangent.
              cNew = cOld;
              pApex = cNew / obj.alpha;

              D(:,:,i) = zeros(6,6);
            end
          end

          % Deviatoric stress is zero at the apex.
          sigma(i,1:6) = (pApex * I)';

          oldCohesion(i,1) = cNew;

        end
      end
    end



    function D = getElasticTensor(obj)
      D = zeros(6);
      D([1 8 15]) = 1 - obj.nu;
      D([2 3 7 9 13 14]) = obj.nu;
      D([22 29 36]) = (1 - 2*obj.nu)/2;
      D = obj.E / ((1 + obj.nu) * (1 - 2*obj.nu)) * D;
    end
  end

  methods (Access = private)
    function readMaterialParameters(obj, varargin)

      default = struct('youngModulus',[], ...
                       'poissonRatio',[], ...
                       'frictionAngle',[], ...
                       'dilatancy',[], ...
                       'cohesion',[], ...
                       'hardeningParameter',0.0);

      params = readInput(default,varargin{:});

      obj.E = params.youngModulus;
      obj.nu = params.poissonRatio;
      obj.phi = params.frictionAngle;
      obj.psi = params.dilatancy;
      obj.co = params.cohesion;
      obj.h = params.hardeningParameter;

      obj.M = obj.nu / (1 - obj.nu);
      obj.cM = (1 + obj.nu) * (1 - 2 * obj.nu) / (obj.E * (1 - obj.nu));

      sinPhi = sin(deg2rad(obj.phi));
      cosPhi = cos(deg2rad(obj.phi));
      sinPsi = sin(deg2rad(obj.psi));

      obj.alpha   = 6.0 * sinPhi / (3.0 - sinPhi);
      obj.beta    = 6.0 * sinPsi / (3.0 - sinPsi);
      obj.epsilon = 6.0 * cosPhi / (3.0 - sinPhi);

    end
  end
end