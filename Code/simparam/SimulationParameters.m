classdef SimulationParameters < handle
  %UNTITLED2 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (SetAccess = private, GetAccess = public)
    itMaxNR = 10
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
    relaxFac = 0
    pTarget
    sTarget = 0.4
    NLSolver = 'Newton'
    verbosity = 2
    goOnBackstep = 0;
    isTimeDependent = true;
  end
  
  methods (Access = public)
    function obj = SimulationParameters(fileName,model)
        [~, ~, extension] = fileparts(fileName);
        switch extension
            case '.xml'
                obj.readXMLFile(fileName,model);
            otherwise
                obj.setSimulationParameters(fileName,model);
        end
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
    function setSimulationParameters(obj,fileName,model)
      fid = fopen(fileName,'r');
      if fid == -1
        error('File %s not opened correctly',fileName);
      end
      %
      % Read the model parameters
      %
      %readSPFParameters(obj,fid,fileName);
      % Single-phase Flow
      if model.isVariabSatFlow()
         readVSFParameters(obj,fid,fileName);
      else
         readSPFParameters(obj,fid,fileName);
      end
      token = SimulationParameters.readToken(fid,fileName);
      if ~strcmp(token,'End')
        error('Missing End statement at the end of Simulation settings file %s',fileName);
      end
      fclose(fid);
    end
    
    function readSPFParameters(obj,fid,fName)
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',false,1,fName)
        obj.tMax = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',false,2,fName)
        obj.dtIni = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',false,3,fName)
        obj.dtMin = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',false,4,fName)
        obj.dtMax = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,5,fName)
        obj.multFac = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,6,fName)
        obj.divFac = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',false,7,fName)
        obj.pTarget = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,8,fName)
        obj.relaxFac = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,9, ...
          fName,[0, 0.5, 1])
        obj.theta = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,10,fName)
        obj.relTol = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,11,fName)
        obj.absTol = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,12, ...
          fName,["2", "Inf"])
        obj.pNorm = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,13,fName)
        obj.itMaxNR = str2double(token);
      end
    end
    
    function readVSFParameters(obj,fid,fName)
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',false,1,fName)
        obj.tMax = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',false,2,fName)
        obj.dtIni = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',false,3,fName)
        obj.dtMin = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',false,4,fName)
        obj.dtMax = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,5,fName)
        obj.multFac = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,6,fName)
        obj.divFac = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',false,7,fName)
        obj.pTarget = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,8,fName)
        obj.relaxFac = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,9, ...
          fName,[0, 0.5, 1])
        obj.theta = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,10,fName)
        obj.relTol = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,11,fName)
        obj.absTol = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,12, ...
          fName,["2", "Inf"])
        obj.pNorm = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'string',true,13, ...
          fName,["Newton", "Picard"])
        obj.NLSolver = token;
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,14,fName)
        obj.itMaxNR = str2double(token);
      end
    end
    
    function readPoromechParameters(obj,fid,fName)
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,1,fName)
        obj.relTol = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,2,fName)
        obj.absTol = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,3, ...
          fName,["2", "Inf"])
        obj.pNorm = str2double(token);
      end
      %
      token = SimulationParameters.readToken(fid,fName);
      if ~SimulationParameters.checkEntry(token,'number',true,4,fName)
        obj.itMaxNR = str2double(token);
      end
      %
      obj.tMax = 1;
      obj.dtIni = 1.e10;
      obj.dtMin = 1;
      obj.dtMax = 1;
    end

    function readXMLFile(obj,fileName,model)
        %READXMLFILE - function to read the simulation parameters file in
        %xml and construct the class object.

        fdata = readstruct(fileName);

        obj.tIni = obj.checkSlotXML(fdata,'Time','Start',0.);
        obj.tMax = obj.checkSlotXML(fdata,'Time','End',1.);
        obj.dtIni = obj.checkSlotXML(fdata,'Time','DtInic',1.);
        obj.dtMin = obj.checkSlotXML(fdata,'Time','DtMin',1.);
        obj.dtMax = obj.checkSlotXML(fdata,'Time','DtMax',1.);

        obj.multFac = obj.checkSlotXML(fdata,'TimeStepControl','IncAugmentation',1.1);
        obj.divFac = obj.checkSlotXML(fdata,'TimeStepControl','IncChopping',2.);
        obj.relaxFac = obj.checkSlotXML(fdata,'TimeStepControl','RelaxationFac',0.);

        obj.pNorm = obj.checkSlotXML(fdata,'Tolerance','PresIncTarget',2);
        obj.absTol = obj.checkSlotXML(fdata,'Tolerance','Absolute',1.e-10);
        obj.relTol = obj.checkSlotXML(fdata,'Tolerance','Relative',1.e-6);
        obj.pTarget = obj.checkSlotXML(fdata,'Tolerance','PresIncTarget',1.e6);

        obj.NLSolver = obj.checkSlotXML(fdata,'Solver','Type','Newton');
        if ~model.isVariabSatFlow()
            obj.NLSolver = 'Newton';
        end
        obj.theta = obj.checkSlotXML(fdata,'Solver','Theta',1.);
        obj.itMaxNR = obj.checkSlotXML(fdata,'Solver','MaxIteration',10);
    end

  end
  
  methods (Static = true)
    function token = readToken(fid,fileName)
      flEof = feof(fid);   % end-of-file flag
      if flEof == 1
        error('End of file while reading Simulation settings file %s',fileName);
      else
        token = strtrim(strtok(fgetl(fid),'%'));
      end  
    end
    
    function isDefault = checkEntry(token,type,defAdm,i,fileName,varargin)
      isDefault = false;
      if isnan(str2double(token))
        if strcmp(token,'Default')
          if defAdm
            isDefault = true;
          else
            error('There is no default value for the parameter in row %d of file %s',i,fileName);
          end
        elseif strcmp(type,'number')
          error('Expected numeric entry in row %d of file %s',i,fileName);
        elseif strcmp(type,'string')
          if ~isempty(varargin)
            if ~ismember(token,varargin{1})
              str = sprintf(' %s',varargin{1});
              fmt = ['Invalid entry in row %d in file %s\nExpected options are:' str];
              error(fmt,i,fileName);
            end
          end
        end
      else
        if strcmp(type,'string')
          error('Expected string entry in row %d in file %s',i,fileName);
        elseif strcmp(type,'number')
          if ~isempty(varargin)
            if ~ismember(token,varargin{1})
              str = sprintf(' %s',varargin{1});
              fmt = ['Invalid entry in row %d in file %s\nExpected options are:' str];
              error(fmt,i,fileName);
            end
          end
        end
      end
    end
  
    function data = checkSlotXML(stname,root,varname,default)
        %CHECKSLOTXML - check the field is defined for the xml struct and
        %return the value or a default value.
        if isfield(stname,root)
            if isfield(stname.(root),varname)
                if strcmp(stname.(root).(varname),"") || strcmp(stname.(root).(varname),"Default")
                    data = default;
                else
                    data = stname.(root).(varname);
                end
            else
                data = default;
            end
        else
            data=default;
        end
    end

  end
end

