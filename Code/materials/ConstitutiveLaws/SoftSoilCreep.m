classdef SoftSoilCreep < handle
  % SOFTSOILCREEP Rate-dependent soft-soil constitutive law (Isotton et al., 2019)
  %
  % Implementation of the strain-driven EPRDM algorithm based on the 
  % Vermeer-Neher soft-soil creep model. Fully compatible with GReS core mechanics.
  %
  % State variables stored at each Gauss point (6-column GReS array):
  %   status(:,1) = reference equivalent pressure p_c,r
  %   status(:,2) = accumulated plastic multiplier / strain gamma
  %   status(:,3:6) = reserved / unused (padded to match GReS core)

  properties (Access = private)
    ni;
    lambda;
    kappa;
    mu;
    tau;
    cohes;
    M;
    fric_ang;
    OCR;
    initialPcR;
    itmax = 50;
    atol = 1.e-12;
    rtol = 1.e-8;
  end

  methods (Access = public)
    % Class constructor method supporting Struct, XML, Key-Value pairs, or legacy handles
    function obj = SoftSoilCreep(varargin)
      obj.readMaterialParameters(varargin{:});
    end

    % Initialize Gauss-point status array (returns 6 columns to match GReS core expectations)
    function [status] = initializeStatus(obj, sigma)
      nptGauss = size(sigma, 1);
      status = zeros(nptGauss, 6);
      for i = 1 : nptGauss
        if ~isempty(obj.initialPcR) && isfinite(obj.initialPcR)
          pc0 = obj.initialPcR;
        else
          % Convert tension-positive GReS stress to positive compression mean pressure
          p0_comp = max(-sum(sigma(i, 1:3)) / 3, 1e-10);
          q0 = sqrt(max(0.5 * ((sigma(i,1)-sigma(i,2))^2 + (sigma(i,2)-sigma(i,3))^2 + (sigma(i,3)-sigma(i,1))^2) + ...
              3 * (sigma(i,4)^2 + sigma(i,5)^2 + sigma(i,6)^2), 0.0));
          
          tanPhi = tan(deg2rad(obj.fric_ang));
          if abs(tanPhi) < 1e-12
            shift = 0.0;
          else
            shift = obj.cohes / tanPhi;
          end
          pHat0 = max(p0_comp + shift, 1e-10);
          pc0 = pHat0 + q0^2 / (obj.M^2 * pHat0);
        end
        status(i, 1) = max(pc0 * obj.OCR, 1e-6); % p_c,r
        status(i, 2) = 0.0;                      % gamma
        status(i, 3:6) = 0.0;                    % padding
      end
    end

    % Material stiffness matrix calculation matching GReS Poromechanics interface
    function [DAll, sigmaOut, statusOut] = getStiffnessMatrix(obj, sigmaIn, depsilon, dt, statusIn, cellID) %#ok<INUSD>
      nptGauss = size(sigmaIn, 1);
      DAll = zeros(6, 6, nptGauss);
      sigmaOut = zeros(size(sigmaIn)); % Preallocated
      
      % Self-healing check for uninitialized or zero status matrices from Poromechanics.initState()
      if isempty(statusIn) || isscalar(statusIn) || size(statusIn, 1) ~= nptGauss || all(statusIn(:) == 0)
        statusClean = obj.initializeStatus(sigmaIn);
      else
        statusClean = statusIn;
        % Self-healing check if p_c,r is uninitialized (0 or negative)
        if any(statusClean(:, 1) <= 0)
          initSt = obj.initializeStatus(sigmaIn);
          statusClean(:, 1) = initSt(:, 1);
        end
      end

      statusOut = zeros(nptGauss, 6); % GReS core expects 6 columns

      for i = 1 : nptGauss
        peqp0_i = statusClean(i, 1);
        gamma0_i = statusClean(i, 2);
        
        [DAll(:,:,i), sigmaOut(i,1:6), peqpOut_i, gammaOut_i] = obj.innerSolverCoord(...
            -sigmaIn(i,1:6)', -depsilon(i,1:6)', dt, peqp0_i, gamma0_i);
            
        statusOut(i, 1) = peqpOut_i;
        statusOut(i, 2) = gammaOut_i;
        statusOut(i, 3:6) = 0.0;
      end
      
      sigmaOut = -sigmaOut;

      for i  = 1:size(DAll,3)
         if norm(DAll(:,:,i)-DAll(:,:,i)','f') > 1e-12
            printf('Error in norm, %d DAll page not symmetric',i);
            error('');
         end
      end
    end

    % Method that returns the M factor
    function m = getMFactor(obj)
      m = obj.M;
    end
  end

  methods (Access = private)
    function readMaterialParameters(obj, varargin)
      default = struct( ...
        'poissonRatio', [], ...
        'modifiedCompressionIndex', [], ...
        'modifiedSwellingIndex', [], ...
        'creepIndex', [], ...
        'referenceTime', [], ...
        'cohesion', 0.0, ...
        'frictionAngle', [], ...
        'K0NC', missing, ...
        'initialOCR', 1.0, ...
        'initialPreconsolidationPressure', missing, ...
        'criticalStateSlope', missing, ...
        'maxIterations', 50, ...
        'absoluteTolerance', 1e-12, ...
        'relativeTolerance', 1e-8);

      if nargin > 2 && isnumeric(varargin{1}) && (ischar(varargin{2}) || isstring(varargin{2}))
        % Legacy positional reading: SoftSoilCreep(fID, matFileName)
        fID = varargin{1};
        matFileName = varargin{2};
        tmpVec = readDataInLine(fID, matFileName, 9);
        obj.ni       = tmpVec(1);
        obj.lambda   = tmpVec(2);
        obj.kappa    = tmpVec(3);
        obj.mu       = tmpVec(4);
        obj.tau      = tmpVec(5);
        obj.cohes    = tmpVec(6);
        obj.M        = tmpVec(7);
        obj.fric_ang = tmpVec(8);
        obj.OCR      = tmpVec(9);
      else
        % Modern XML / Struct / Key-Value parsing
        params = readInput(default, varargin{:});
        
        obj.ni       = params.poissonRatio;
        obj.lambda   = params.modifiedCompressionIndex;
        obj.kappa    = params.modifiedSwellingIndex;
        obj.mu       = params.creepIndex;
        obj.tau      = params.referenceTime;
        obj.cohes    = params.cohesion;
        obj.fric_ang = params.frictionAngle;
        obj.OCR      = params.initialOCR;
        
        if isfield(params, 'initialPreconsolidationPressure') && ...
           ~isempty(params.initialPreconsolidationPressure) && ...
           ~ismissing(params.initialPreconsolidationPressure)
          obj.initialPcR = params.initialPreconsolidationPressure;
        end

        if isfield(params, 'maxIterations'), obj.itmax = params.maxIterations; end
        if isfield(params, 'absoluteTolerance'), obj.atol = params.absoluteTolerance; end
        if isfield(params, 'relativeTolerance'), obj.rtol = params.relativeTolerance; end

        if ~ismissing(params.criticalStateSlope) && ~isempty(params.criticalStateSlope)
          obj.M = params.criticalStateSlope;
        elseif ~isempty(params.K0NC) && ~ismissing(params.K0NC)
          r = obj.lambda / obj.kappa;
          k0 = params.K0NC;
          denom = (1+2*k0)*(1-2*obj.ni)*r - (1-k0)*(1+obj.ni);
          obj.M = 3 * sqrt(((1-k0)^2)/(1+2*k0)^2 + ((1-k0)*(1-2*obj.ni)*(r-1))/denom);
        else
          error('SoftSoilCreep:MissingCriticalStateSlope', ...
            'Either criticalStateSlope or K0NC must be provided.');
        end
      end
    end

    function [DAll, sigmaOut, peqpOut, gammaOut] = innerSolverCoord(obj, sigmaIn, deps, dt, peqp0In, gamma0In)
      nTrials = 4;
      nSubs = 1;
      for iTrial = 1 : nTrials
        depsLoc = deps / nSubs;
        dtLoc = dt / nSubs;
        peqpLoc = peqp0In;
        gammaLoc = gamma0In;
        sigmaInLoc = sigmaIn;
        iSub = 1;
        flag = true;
        while (iSub <= nSubs && flag)
          [p, q, peqpOut, sigmaUpd, DAll, flag] = ...
              solve2(sigmaInLoc, depsLoc, peqpLoc, obj.ni, obj.lambda, obj.kappa, obj.mu, obj.tau, ...
              obj.cohes, obj.M, obj.fric_ang, dtLoc, obj.itmax, obj.atol, obj.rtol);
          
          % Accumulate plastic multiplier gamma using the rate-dependent creep formula
          tanPhi = tan(deg2rad(obj.fric_ang));
          shift = obj.cohes / max(tanPhi, 1e-12);
          phat = max(p + shift, 1e-10);
          peq = phat + q^2 / (obj.M^2 * phat);
          dpeq_dp = 1.0 - q^2 / (obj.M^2 * phat^2);
          beta = sign(dpeq_dp) * max(abs(dpeq_dp), 1e-4);
          
          ratio = min(max(peq / max(peqpOut, 1e-12), 1e-12), 1e10);
          exponent = (obj.lambda - obj.kappa) / max(obj.mu, 1e-12);
          creepRate = (obj.mu / obj.tau) * (1.0 / beta) * (ratio^exponent);
          dGamma = dtLoc * creepRate;
          
          sigmaInLoc = sigmaUpd;
          peqpLoc = peqpOut;
          gammaLoc = gammaLoc + dGamma;
          iSub = iSub + 1;
        end
        if (flag)
          break;
        end
        nSubs = nSubs * 2;
      end
      
      if ~flag
        % Elastic fallback stiffness matrix protecting global Newton iteration
        p_comp = max((sigmaIn(1)+sigmaIn(2)+sigmaIn(3))/3, 1e-10);
        E = 3 * (1 - 2*obj.ni) * p_comp / obj.kappa;
        DAll = zeros(6);
        DAll([1 8 15]) = 1 - obj.ni;
        DAll([2 3 7 9 13 14]) = obj.ni;
        DAll([22 29 36]) = (1 - 2*obj.ni) / 2;
        DAll = (E / ((1 + obj.ni) * (1 - 2*obj.ni))) * DAll;
        sigmaOut = sigmaIn + DAll * deps;
        peqpOut = peqp0In;
        gammaOut = gamma0In;
      else
        sigmaOut = sigmaUpd;
        gammaOut = gammaLoc;
      end
    end
  end
end