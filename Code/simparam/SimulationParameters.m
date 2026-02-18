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
 
    function readSimulationParameters(obj,fileName)
        %READXMLFILE - function to read the simulation parameters file in
        %xml and construct the class object.

        input = readstruct(fileName,AttributeSuffix="");

        if isfield(input,'simParam')
          input = input.simParam;
        end

        if isfield(input,"fileName")
          assert(isscalar(fieldnames(input)),"FileName, " + ...
            " must be a unique parameter.");
          input = readstruct(input.fileName,AttributeSuffix="");
        end

        time = input.("Time");

        obj.tIni = getXMLData(time,0,'Start');
        obj.tMax = getXMLData(time,[],'End');
        obj.dtIni = getXMLData(time,[],'DtInit');
        obj.dtMin = getXMLData(time,[],'DtMin');
        obj.dtMax = getXMLData(time,[],'DtMax');
        obj.multFac = getXMLData(time,1.1,'incrementFactor');
        obj.divFac = getXMLData(time,2.,'choppingFactor');

        if isfield(input,"Solver")
          solver = input.("Solver");
          obj.absTol = getXMLData(solver,1e-10,'AbsoluteTolerance');
          obj.relTol = getXMLData(solver,1e-6','RelativeTolerance');
          %obj.theta = getXMLData(solver,1.,'Theta');
          obj.itMaxNR = getXMLData(solver,10,'MaxNLIteration');
          obj.itMaxConfig = getXMLData(solver,10,'MaxConfigurationIteration');
          obj.attemptSimplestConfiguration = ...
            logical(getXMLData(solver,0,'resetConfiguration'));
        end
    end

  end

end

