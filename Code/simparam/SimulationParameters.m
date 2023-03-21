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
    tMax
    multFac = 1.1
    divFac = 2
    relaxFac = 0
    pTarget
  end
  
  methods (Access = public)
    function obj = SimulationParameters(fileName)
      obj.setSimulationParameters(fileName);
    end
  end
  
  methods (Access = private)
    function setSimulationParameters(obj,fileName)
      fid = fopen(fileName,'r');
      if fid == -1
        error('File %s not opened correctly',fileName);
      end
      block = '';
      [flEof,line] = SimulationParameters.readLine(fid);
      while (~strcmpi(line,'end'))
        if flEof == 1
          error('End of file while reading Simulation settings file %s',obj.fileName);
        end
        line = strtrim(strtok(line,'%'));
        block = [block, line, ' '];
        [flEof,line] = SimulationParameters.readLine(fid);
      end
      fclose(fid);
      blockSplt = strsplit(string(deblank(block)));
      tmpVec = str2double(blockSplt);
      NaNParam = isnan(tmpVec);
      if ~all(strcmpi(blockSplt(NaNParam),'default'))
        error('There are invalid entries in %s file',fileName);
      end
      %
      if NaNParam(1) == 0
        obj.tMax = tmpVec(1);
      else
        error('There is no default value for tMax - row 1 of %s file',fileName);
      end
      %
      if NaNParam(2) == 0
        obj.dtIni = tmpVec(2);
      else
        error('There is no default value for dtIni - row 2 of %s file',fileName);
      end
      %
      if NaNParam(3) == 0
        obj.dtMin = tmpVec(3);
      else
        error('There is no default value for dtMin - row 3 of %s file',fileName);
      end
      %
      if NaNParam(4) == 0
        obj.dtMax = tmpVec(4);
      else
        error('There is no default value for dtMax - row 4 of %s file',fileName);
      end
      if NaNParam(5) == 0; obj.multFac = tmpVec(5); end
      if NaNParam(6) == 0; obj.divFac = tmpVec(6); end
      %
      if NaNParam(7) == 0
        obj.pTarget = tmpVec(7);
      else
        error('There is no default value for pTarget - row 7 of %s file',fileName);
      end
      if NaNParam(8) == 0; obj.relaxFac = tmpVec(8); end
      if NaNParam(9) == 0; obj.theta = tmpVec(9); end
      if NaNParam(10) == 0; obj.relTol = tmpVec(10); end
      if NaNParam(11) == 0; obj.absTol = tmpVec(11); end
      %
      if NaNParam(12) == 0
        if ismember(tmpVec(12),[2 Inf])
          obj.pNorm = tmpVec(12);
        else
          error('Invalid value for pTarget - row 7 of %s file\nAccepted p-norms are: 2 and Inf',fileName);
        end
      end
      if NaNParam(13) == 0; obj.itMaxNR = tmpVec(13); end
    end
  end
  
  methods (Static = true)
    function [flEof,line] = readLine(fid)
      flEof = feof(fid);   % end-of-file flag
      if flEof == 1
        line = '';
      else
        line = deblank(fgetl(fid));
        if isempty(line)
          error('No blank lines are admitted in Simulatio parameters file');
        end
      end  
    end
  end
end

