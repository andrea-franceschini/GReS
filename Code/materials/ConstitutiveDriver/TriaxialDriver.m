classdef TriaxialDriver < handle
  % General class to run simple triaxial tests using a specified
  % constitutive law.
  %
  % Voigt convention used by the driver:
  %   stress/strain = [0, 1, 2, 12, 02, 01]
  % where direction 0 is axial and directions 1 and 2 are radial.
  % Therefore, triaxial symmetry means (in MATLAB indexing)
  %   eps(3) = eps(2), sig(3) = sig(2).
  % Inspired by the GEOS triaxial driver

  properties
    constLaw          % Constitutive law object
    func = struct()   % Prescribed loading functions
    params            % Driver parameters
    control           % "StressControl", "StrainControl", or "MixedControl"
    nStep             % Number of loading steps
    outFile           % Optional output file
  end

  properties (Access = private)
    % Linear-index order is intentional:
    %   funFlag(1) = axialStress  -> funFlag(1,1)
    %   funFlag(2) = radialStress -> funFlag(2,1)
    %   funFlag(3) = axialStrain  -> funFlag(1,2)
    %   funFlag(4) = radialStrain -> funFlag(2,2)
    % Equivalently, rows are axial/radial and columns are stress/strain.
    funFlag = false(2,2)
  end

  properties (SetAccess = private, GetAccess = public)
    % index for output table
    TIME = 1
    SIG0 = 2
    SIG1 = 3
    SIG2 = 4
    EPS0 = 5
    EPS1 = 6
    EPS2 = 7
    ITER = 8
    NORM = 9
    EXITFLAG = 10
    NUMCOLS = 10
  end

  methods
    function obj = TriaxialDriver(varargin)

      default = struct('nStep',1,'constLaw',[],'outFile',missing);
      parm = readInput(default,varargin{:});

      obj.nStep = parm.nStep;
      obj.constLaw = parm.constLaw;
      obj.outFile = parm.outFile;

      % default parameters
      obj.params = struct();
      obj.params.initialStress = 0.0;
      obj.params.maxIter       = 25;
      obj.params.maxCuts       = 8;
      obj.params.newtonTol     = 1.0e-10;
      obj.params.minJacobian   = 1.0e-30;

    end

    function setFunction(obj, type, xval, yval)
      types = ["axialStress", "radialStress", "axialStrain", "radialStrain"];
      typeID = strcmpi(type, types);

      if ~any(typeID)
        error(['Unrecognized function type. Available types are ', ...
               '"axialStress", "radialStress", "axialStrain", "radialStrain".'])
      end

      typeID = find(typeID, 1);

      if obj.funFlag(typeID)
        warning("Function %s already defined. It will be overwritten.", types(typeID));
      end

      xval = xval(:);
      yval = yval(:);

      if numel(xval) ~= numel(yval)
        error('The x and y arrays for %s must have the same number of entries.', types(typeID));
      end
      if numel(xval) < 2
        error('The function %s must contain at least two tabulated points.', types(typeID));
      end
      if any(diff(xval) <= 0)
        error('The x values for %s must be strictly increasing.', types(typeID));
      end

      obj.funFlag(typeID) = true;
      obj.func.(types(typeID)).x = xval;
      obj.func.(types(typeID)).y = yval;
    end

    function setParameters(obj, varargin)

      % Override default parameters
      obj.params = readInput(obj.params,varargin{:});

    end

    function out = launch(obj) 
      % Return output in matrix format:
      % [time, sig0, sig1, sig2, eps0, eps1, eps2, iter, norm, exitflag]

      obj.validateControlFlags();
      obj.control = obj.detectControlMode();

      % prepare output table
      out = obj.initializeTable(obj.control);

      c = obj.constLaw.initializeStatus(0);


      switch obj.control
        case "StressControl"
          out = obj.runStressControl(out,c);
        case "StrainControl"
          out = obj.runStrainControl(out,c);
        case "MixedControl"
          out = obj.runMixedControl(out,c);
        otherwise
          error('Unknown control mode.');
      end

      if ~ismissing(obj.outFile)
        writematrix(out, obj.outFile);
      end
    end
  end

  methods (Access = private)


    function out = runStrainControl(obj, out, status)

      sigma = zeros(6,1);
      sigma(1:3) = obj.params.initialStress;

      for n = 2:size(out,1)
        strainIncrement = zeros(6,1);
        strainIncrement(1) = out(n, obj.EPS0) - out(n-1, obj.EPS0);
        strainIncrement(2) = out(n, obj.EPS1) - out(n-1, obj.EPS1);
        strainIncrement(3) = out(n, obj.EPS2) - out(n-1, obj.EPS2);

        % timeIncrement = out(n, obj.TIME) - out(n-1, obj.TIME);
        [~, sigma, status] = obj.localConstitutiveUpdate(sigma, strainIncrement, status);

        out(n, obj.SIG0) = sigma(1);
        out(n, obj.SIG1) = sigma(2);
        out(n, obj.SIG2) = sigma(3);
        out(n, obj.ITER) = 0;
        out(n, obj.NORM) = 0;
        out(n, obj.EXITFLAG) = 1;

      end
    end

    function out = runStressControl(obj, out, oldStatus)

      sigma = zeros(6,1);
      sigma(1:3) = obj.params.initialStress;

      maxIter = obj.params.maxIter;
      maxCuts = obj.params.maxCuts;
      tol     = obj.params.newtonTol;
      scale   = obj.computeStressScale(out, [obj.SIG0, obj.SIG1, obj.SIG2]);

      for n = 2:size(out,1)
        %timeIncrement = out(n, obj.TIME) - out(n-1, obj.TIME);
        strainIncrement = zeros(6,1);
        deltaStrainIncrement = zeros(2,1);

        target = [out(n, obj.SIG0); out(n, obj.SIG1)];
        normR = inf;
        normZero = inf;
        cuts = 0;
        exitflag = 0;

        for k = 0:maxIter-1
          [stiffness, trialStress, status] = obj.localConstitutiveUpdate(sigma, strainIncrement, oldStatus);

          resid = scale .* ([trialStress(1); trialStress(2)] - target);
          normR = norm(resid);

          if k == 0
            normZero = normR;
          end

          if normR < tol
            sigma = trialStress;
            exitflag = 1;
            break
          elseif k > 0 && normR > normZero && cuts < maxCuts
            cuts = cuts + 1;
            deltaStrainIncrement = 0.5 .* deltaStrainIncrement;
            strainIncrement(1) = strainIncrement(1) + deltaStrainIncrement(1);
            strainIncrement(2) = strainIncrement(2) + deltaStrainIncrement(2);
            strainIncrement(3) = strainIncrement(2);
          else
            cuts = 0;
            jacobian = scale .* [stiffness(1,1), stiffness(1,2) + stiffness(1,3); ...
                                 stiffness(2,1), stiffness(2,2) + stiffness(2,3)];

            if rcond(jacobian) < obj.params.minJacobian
              exitflag = -2;
              sigma = trialStress;
              break
            end

            deltaStrainIncrement = jacobian \ resid;

            strainIncrement(1) = strainIncrement(1) - deltaStrainIncrement(1);
            strainIncrement(2) = strainIncrement(2) - deltaStrainIncrement(2);
            strainIncrement(3) = strainIncrement(2);
          end
        end

        oldStatus = status;

        out(n, obj.EPS0) = out(n-1, obj.EPS0) + strainIncrement(1);
        out(n, obj.EPS1) = out(n-1, obj.EPS1) + strainIncrement(2);
        out(n, obj.EPS2) = out(n, obj.EPS1);
        out(n, obj.ITER) = k;
        out(n, obj.NORM) = normR;
        out(n, obj.EXITFLAG) = exitflag;

        if exitflag ~= 1
          break
        end
      end
    end

    function out = runMixedControl(obj, out, oldStatus)

      sigma = zeros(6,1);
      sigma(1:3) = obj.params.initialStress;

      maxIter = obj.params.maxIter;
      maxCuts = obj.params.maxCuts;
      tol     = obj.params.newtonTol;
      scale   = obj.computeStressScale(out, [obj.SIG1, obj.SIG2]);

      for n = 2:size(out,1)
        %timeIncrement = out(n, obj.TIME) - out(n-1, obj.TIME);
        strainIncrement = zeros(6,1);
        strainIncrement(1) = out(n, obj.EPS0) - out(n-1, obj.EPS0);

        deltaRadialStrain = 0;
        targetRadialStress = out(n, obj.SIG1);
        normR = inf;
        normZero = inf;
        cuts = 0;
        exitflag = 0;

        for k = 0:maxIter-1
          [stiffness, trialStress, status] = obj.localConstitutiveUpdate(sigma, strainIncrement, oldStatus);

          resid = scale .* (trialStress(2) - targetRadialStress);
          normR = abs(resid);

          if k == 0
            normZero = normR;
          end

          if normR < tol
            sigma = trialStress;
            exitflag = 1;
            break
          elseif k > 0 && normR > normZero && cuts < maxCuts
            cuts = cuts + 1;
            deltaRadialStrain = 0.5 * deltaRadialStrain;
            strainIncrement(2) = strainIncrement(2) + deltaRadialStrain;
            strainIncrement(3) = strainIncrement(2);
          else
            cuts = 0;
            jacobian = scale .* (stiffness(2,2) + stiffness(2,3));

            if abs(jacobian) < obj.params.minJacobian
              exitflag = -2;
              sigma = trialStress;
              break
            end

            deltaRadialStrain = jacobian \ resid;
            strainIncrement(2) = strainIncrement(2) - deltaRadialStrain;
            strainIncrement(3) = strainIncrement(2);
          end
        end

        % commit hardening variable
        oldStatus = status;

        out(n, obj.SIG0) = sigma(1);
        out(n, obj.SIG1) = sigma(2);
        out(n, obj.SIG2) = sigma(3);
        out(n, obj.EPS1) = out(n-1, obj.EPS1) + strainIncrement(2);
        out(n, obj.EPS2) = out(n, obj.EPS1);
        out(n, obj.ITER) = k;
        out(n, obj.NORM) = normR;
        out(n, obj.EXITFLAG) = exitflag;

        if exitflag ~= 1
          break
        end
      end
    end

    function out = initializeTable(obj, mode)

      axialName = obj.getAxialFunctionName();
      axialFun = obj.func.(axialName);

      minTime = axialFun.x(1);
      maxTime = axialFun.x(end);
      time = linspace(minTime, maxTime, obj.nStep + 1).';

      out = zeros(obj.nStep + 1, obj.NUMCOLS);
      out(:, obj.TIME) = time;
      out(1, obj.SIG0:obj.SIG2) = obj.params.initialStress;
      out(:, obj.EXITFLAG) = NaN;
      out(:, obj.ITER) = NaN;
      out(:, obj.NORM) = NaN;

      switch mode
        case "StressControl"
          fAxi = obj.func.axialStress;
          fRad = obj.func.radialStress;
          axi = interp1(fAxi.x, fAxi.y, time, 'linear', 'extrap');
          rad = interp1(fRad.x, fRad.y, time, 'linear', 'extrap');
          out(:, obj.SIG0) = axi;
          out(:, obj.SIG1) = rad;
          out(:, obj.SIG2) = rad;

        case "StrainControl"
          fAxi = obj.func.axialStrain;
          fRad = obj.func.radialStrain;
          axi = interp1(fAxi.x, fAxi.y, time, 'linear', 'extrap');
          rad = interp1(fRad.x, fRad.y, time, 'linear', 'extrap');
          out(:, obj.EPS0) = axi;
          out(:, obj.EPS1) = rad;
          out(:, obj.EPS2) = rad;

        case "MixedControl"
          fAxi = obj.func.axialStrain;
          fRad = obj.func.radialStress;
          axi = interp1(fAxi.x, fAxi.y, time, 'linear', 'extrap');
          rad = interp1(fRad.x, fRad.y, time, 'linear', 'extrap');
          out(:, obj.EPS0) = axi;
          out(:, obj.SIG1) = rad;
          out(:, obj.SIG2) = rad;
      end

      obj.checkInitialStressConsistency(out, mode);
      out(1, obj.EXITFLAG) = 1;
      out(1, obj.ITER) = 0;
      out(1, obj.NORM) = 0;
    end

    function [stiffness, sigmaOut, status] = localConstitutiveUpdate(obj, sigma, strainIncrement, oldStatus)
      if isempty(obj.constLaw)
        error('No constitutive law object was provided.');
      end

      % Expected user-side interface, following the commented skeleton:
      %   [constMat, sigmaOut, ...] = constLaw.getStiffnessMatrix(sigma, strainIncrement, 0, 0, timeIncrement)
      [stiffness, sigmaOut, status] = obj.constLaw.getStiffnessMatrix(sigma', strainIncrement', 0, oldStatus, 1);

      sigmaOut = sigmaOut(:);
      if numel(sigmaOut) ~= 6
        error('The constitutive law must return a 6-component stress vector.');
      end
      if ~isequal(size(stiffness), [6, 6])
        error('The constitutive law must return a 6-by-6 stiffness matrix.');
      end
    end

    function scale = computeStressScale(obj, out, cols)
      vals = out(2:end, cols);
      vals = vals(isfinite(vals));
      denom = mean(abs(vals));
      if isempty(denom) || denom <= eps
        scale = 1.0;
      else
        scale = 1.0 / denom;
      end
    end

    function validateControlFlags(obj)
      if any(sum(obj.funFlag, 2) == 2)
        error('Cannot set both stress and strain in the same direction.');
      end
      if any(sum(obj.funFlag, 2) == 0)
        error('Missing function for axial or radial direction.');
      end
    end

    function mode = detectControlMode(obj)
      if obj.funFlag(1,1) && obj.funFlag(2,1)
        mode = "StressControl";
      elseif obj.funFlag(1,2) && obj.funFlag(2,2)
        mode = "StrainControl";
      elseif obj.funFlag(1,2) && obj.funFlag(2,1)
        mode = "MixedControl";
      else
        error(['Unsupported mixed mode. This driver supports only ', ...
               'axialStrain + radialStress for MixedControl.']);
      end
    end

    function name = getAxialFunctionName(obj)
      if obj.funFlag(1,1)
        name = "axialStress";
      elseif obj.funFlag(1,2)
        name = "axialStrain";
      else
        error('Missing axial function.');
      end
    end

    function checkInitialStressConsistency(obj, out, mode)
      if mode == "StressControl"
        if abs(obj.params.initialStress - out(1, obj.SIG0)) > 1.0e-6
          error('Initial stress is inconsistent with axialStress at the initial time.');
        end
        if abs(obj.params.initialStress - out(1, obj.SIG1)) > 1.0e-6
          error('Initial stress is inconsistent with radialStress at the initial time.');
        end
      elseif mode == "MixedControl"
        if abs(obj.params.initialStress - out(1, obj.SIG1)) > 1.0e-6
          error('Initial stress is inconsistent with radialStress at the initial time.');
        end
      end
    end

  end
end
