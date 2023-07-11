classdef SSCM < handle
  % ELASTIC ISOTROPIC material class

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
    itmax = 50;
    atol = 1.e-12;
    rtol = 1.e-8;
  end

  methods (Access = public)
    % Class constructor method
    function obj = SSCM(fID, matFileName)
      % Calling the function to set the object properties
      obj.readMaterialParameters(fID, matFileName);
    end

    % Initialize status
    function [status] = initializeStatus(obj, sigma)
      nptGauss = size(sigma,1);
      status = zeros(nptGauss,2);
      for i = 1 : nptGauss
        p0 = (sigma(i,1)+sigma(i,2)+sigma(i,3))/3;
        q0 = sqrt(sigma(i,1)*(sigma(i,1)-sigma(i,2))+sigma(i,2)*(sigma(i,2)-sigma(i,3))+ ...
            sigma(i,3)*(sigma(i,3)-sigma(i,1))+3*(sigma(i,4)^2+sigma(i,5)^2+sigma(i,6)^2));
        status(i,1) = -(p0 + q0^2/(obj.M^2*p0))*obj.OCR;
      end
    end
    %
    % Material stiffness matrix calculation using the object properties
    function [DAll, sigmaOut, peqpOut] = getStiffnessMatrix(obj, sigmaIn, depsilon, dt, peqp0)
      nptGauss = size(sigmaIn,1);
      DAll = zeros(6,6,nptGauss);
      peqpOut = peqp0;
      for i = 1 : nptGauss
        [DAll(:,:,i), sigmaOut(i,:), peqpOut(i,1)] = innerSolver(obj, -sigmaIn(i,:)', -depsilon(i,:)', dt, peqp0(i,1));
      end
      sigmaOut = -sigmaOut;
    end

    % Method that returns the M factor
    function m = getMFactor(obj)
      m = obj.M;
    end

    % Get vertical compressibility
%     function cM = getRockCompressibility(obj)
%       cM = obj.cM;
%     end
  end

  methods (Access = private)
    % Assigning material parameters (check also the Materials class)
    % to object properties
    function readMaterialParameters(obj, fID, matFileName)
      tmpVec = readDataInLine(fID, matFileName, 9);
      obj.ni = tmpVec(1);
      obj.lambda = tmpVec(2);
      obj.kappa = tmpVec(3);
      obj.mu = tmpVec(4);
      obj.tau = tmpVec(5);
      obj.cohes = tmpVec(6);
      obj.M = tmpVec(7);
      obj.fric_ang = tmpVec(8);
      obj.OCR = tmpVec(9);
    end

    function [DAll, sigmaOut, peqpOut] = innerSolver(obj, sigmaIn, deps, dt, peqp0In)
      nTrials = 2;
      nSubs = 1;
      for iTrial = 1 : nTrials
        depsLoc = deps/nSubs;
        dtLoc = dt/nSubs;
        peqp0 = peqp0In;
        sigmaInLoc = sigmaIn;
        iSub = 1;
        flag = true;
        while (iSub <= nSubs && flag)
          [p, q, peqpOut, sigmaUpd, DAll, flag] = ...
              solve2(sigmaInLoc, depsLoc, peqp0In, obj.ni, obj.lambda, obj.kappa, obj.mu, obj.tau, ...
              obj.cohes, obj.M, obj.fric_ang, dtLoc, obj.itmax, obj.atol, obj.rtol);
          sigmaInLoc = sigmaUpd;
          peqpIn = peqpOut;
          iSub = iSub + 1;
        end
        if (flag)
          break;
        end
        nSubs = nSubs * 2;
      end
      sigmaOut = sigmaUpd;
    end
  end
end
