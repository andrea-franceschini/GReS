classdef SimulationParameters < handle
  %UNTITLED2 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (SetAccess = private, GetAccess = public)
    itMaxNR = 10
    itMaxConfig = 10
    relTol = 1.e-6
    absTol = 1.e-10
    pNorm = 2
    theta = 1
    dtIni
    dtMin
    dtMax
    tIni = 0.
    tMax
    multFac = 1.1
    divFac = 2
    pTarget
    sTarget = 0.4
    NLSolver = 'Newton'
    goOnBackstep = 0;
    isTimeDependent = true;
    attemptSimplestConfiguration = false;
  end

  methods (Access = public)
    function obj = SimulationParameters(fileName,varargin)
      obj.readSimulationParameters(fileName);
    end

    function status = isNewtonNLSolver(obj)
      status = false;
      if strcmp(obj.NLSolver,'Newton')
        status = true;
      end
    end
    
    function status = isPicardNLSolver(obj)
      status = false;
      if strcmp(obj.NLSolver,'Picard')
        status = true;
      end
    end

    function setTimeDependence(obj,flag)
       obj.isTimeDependent = flag;
    end

  end
  
  methods (Access = private)
 
    function readSimulationParameters(obj,varargin)
      %READXMLFILE - function to read the simulation parameters file in
      %xml and construct the class object.

      default = struct( ...
        ... % time params
        "Start", 0.0, ...
        "End", [], ...
        "DtInit", [], ...
        "DtMin", [], ...
        "DtMax", [], ...
        "incrementFactor", 1.1, ...
        "choppingFactor", 2.0, ...
        ...% solver params
        "AbsoluteTolerance", 1e-10, ...
        "RelativeTolerance", 1e-6, ...
        "MaxNLIteration", 10, ...
        "MaxConfigurationIteration", 10, ...
        "resetConfiguration", 0 ...
        );

      params = readInput(default,varargin{:});

      obj.tIni   = params.Start;
      obj.tMax   = params.End;
      obj.dtIni  = params.DtInit;
      obj.dtMin  = params.DtMin;
      obj.dtMax  = params.DtMax;
      obj.multFac = params.incrementFactor;
      obj.divFac  = params.choppingFactor;

      obj.absTol     = params.AbsoluteTolerance;
      obj.relTol     = params.RelativeTolerance;
      obj.itMaxNR    = params.MaxNLIteration;
      obj.itMaxConfig = params.MaxConfigurationIteration;
      obj.attemptSimplestConfiguration = logical(params.resetConfiguration);
    end

  end

end

