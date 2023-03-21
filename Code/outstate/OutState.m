classdef OutState < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = public)
    timeList
  end
  
  properties (Access = private)
    model
    printPath
    mesh
    timeID = 1
    VTK
%     flPrint = true
  end
  
  methods (Access = public)
    function obj = OutState(symMod,mesh,fileName,varargin)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      nIn = nargin;
      data = varargin;
      obj.setOutState(nIn,symMod,mesh,fileName,data)
    end
    
    function printState(obj,stateOld,stateNew)
      if nargin == 2
        if isCoupFlowPoro(obj.model)
          [avStressOld,avStrainOld] = finalizeStatePoro(stateOld);
          [fluidPotOld] = finalizeStateFlow(stateOld);
          printProp = struct('time',stateOld.t,'displ',stateOld.displ, ...
            'stress',avStressOld,'strain',avStrainOld, ...
            'pressure',stateOld.pressure,'potential',fluidPotOld);
        elseif isPoromechanics(obj.model)
          [avStressOld,avStrainOld] = finalizeStatePoro(stateOld);
          printProp = struct('time',stateOld.t,'displ',stateOld.displ, ...
            'stress',avStressOld,'strain',avStrainOld);
        elseif isSinglePhaseFlow(obj.model)
          [fluidPotOld] = finalizeStateFlow(stateOld);
          printProp = struct('time',stateOld.t,'pressure',stateOld.pressure, ...
          'potential',fluidPotOld);
        end
        buildPrintStruct(obj,printProp);
      elseif nargin == 3
        if obj.timeID <= length(obj.timeList)
          while (obj.timeList(obj.timeID) <= stateNew.t)
            assert(obj.timeList(obj.timeID) > stateOld.t, ...
            'Print time %f out of range (%f - %f)',obj.timeList(obj.timeID), ...
            stateOld.t,stateNew.t);
            assert(stateNew.t - stateOld.t > eps('double'),'Dt too small for printing purposes');
            %
            % Linear interpolation
            fac = (obj.timeList(obj.timeID) - stateOld.t)/(stateNew.t - stateOld.t);
            if isCoupFlowPoro(obj.model)
              [avStressOld,avStrainOld] = finalizeStatePoro(stateOld);
              [avStressNew,avStrainNew] = finalizeStatePoro(stateNew);
              [fluidPotOld] = finalizeStateFlow(stateOld);
              [fluidPotNew] = finalizeStateFlow(stateNew);
              printProp = struct('time',obj.timeList(obj.timeID), ...
                'displ',stateNew.displ*fac+stateOld.displ*(1-fac), ...
                'stress',avStressNew*fac+avStressOld*(1-fac), ...
                'strain',avStrainNew*fac+avStrainOld*(1-fac), ...
                'pressure',stateNew.pressure*fac+stateOld.pressure*(1-fac), ...
                'potential',fluidPotNew*fac+fluidPotOld*(1-fac));
            elseif isPoromechanics(obj.model)
              [avStressOld,avStrainOld] = finalizeStatePoro(stateOld);
              [avStressNew,avStrainNew] = finalizeStatePoro(stateNew);
              printProp = struct('time',obj.timeList(obj.timeID), ...
                'displ',stateNew.displ*fac+stateOld.displ*(1-fac), ...
                'stress',avStressNew*fac+avStressOld*(1-fac), ...
                'strain',avStrainNew*fac+avStrainOld*(1-fac));
            elseif isSinglePhaseFlow(obj.model)
              [fluidPotOld] = finalizeStateFlow(stateOld);
              [fluidPotNew] = finalizeStateFlow(stateNew);
              printProp = struct('time',obj.timeList(obj.timeID), ...
                'pressure',stateNew.pressure*fac+stateOld.pressure*(1-fac), ...
                'potential',fluidPotNew*fac+fluidPotOld*(1-fac));
            end
            buildPrintStruct(obj,printProp);
            obj.timeID = obj.timeID + 1;
            if obj.timeID > length(obj.timeList)
              break
            end
          end
        end
      end
    end
  end
  
  methods (Access = private)
    function setOutState(obj,nIn,symMod,mesh,fileName,data)
      if nIn > 3
        obj.printPath = data{1};
      end
      obj.model = symMod;
      obj.mesh = mesh;
      %
      fid = fopen(fileName,'r');
      [flEof,line] = OutState.readLine(fid);
      %
      block = '';
      while ~strcmpi(line,'End')
        line = strtrim(line);
        if isempty(line)
          error('Blank line encountered while reading the print times in file %s',fileName);
        elseif flEof == 1
          error('End of file while reading the print times in file %s',fileName);
        end
        block = [block, line, ' '];
        [flEof,line] = OutState.readLine(fid);
      end
      blockSplt = strsplit(string(deblank(block)));
      if ~(blockSplt == "")
        obj.timeList = str2double(blockSplt);
        if any(isnan(obj.timeList))
          error('There are invalid entries in the list of output times')
        end
      else
        obj.timeList = [0 1];
      end
      fclose(fid);
      %
      obj.VTK = VTKOutput(obj.mesh);
    end
    %
    
    function buildPrintStruct(obj,printProp)
      if isCoupFlowPoro(obj.model)
        nPointProp = 5;
        nCellProp = 12;
      elseif isPoromechanics(obj.model)
        nPointProp = 3;
        nCellProp = 12;
      elseif isSinglePhaseFlow(obj.model)
        nPointProp = 2;
      end
      %
      pointData3D = repmat(struct('name', 1, 'data', 1), nPointProp, 1);
      if isPoromechanics(obj.model)
        cellData3D = repmat(struct('name', 1, 'data', 1), nCellProp, 1);
      end
      %
      if isPoromechanics(obj.model)
        %
        % Displacement
        pointData3D(1).name = 'ux';
        pointData3D(1).data = printProp.displ(1:3:length(printProp.displ));
        pointData3D(2).name = 'uy';
        pointData3D(2).data = printProp.displ(2:3:length(printProp.displ));
        pointData3D(3).name = 'uz';
        pointData3D(3).data = printProp.displ(3:3:length(printProp.displ));
        %
        % Stress
        cellData3D(1).name = 'sx';
        cellData3D(1).data = printProp.stress(:,1);
        cellData3D(2).name = 'sy';
        cellData3D(2).data = printProp.stress(:,2);
        cellData3D(3).name = 'sz';
        cellData3D(3).data = printProp.stress(:,3);
        cellData3D(4).name = 'txy';
        cellData3D(4).data = printProp.stress(:,4);
        cellData3D(5).name = 'tyz';
        cellData3D(5).data = printProp.stress(:,5);
        cellData3D(6).name = 'txz';
        cellData3D(6).data = printProp.stress(:,6);
        %
        % Strain
        cellData3D(7).name = 'ex';
        cellData3D(7).data = printProp.strain(:,1);
        cellData3D(8).name = 'ey';
        cellData3D(8).data = printProp.strain(:,2);
        cellData3D(9).name = 'ez';
        cellData3D(9).data = printProp.strain(:,3);
        cellData3D(10).name = 'gxy';
        cellData3D(10).data = printProp.strain(:,4);
        cellData3D(11).name = 'gyz';
        cellData3D(11).data = printProp.strain(:,5);
        cellData3D(12).name = 'gxz';
        cellData3D(12).data = printProp.strain(:,6);
      end
      %
      if isSinglePhaseFlow(obj.model)
        %
        % Pressure
        pointData3D(nPointProp-1).name = 'press';
        pointData3D(nPointProp-1).data = printProp.pressure;
        pointData3D(nPointProp).name = 'potential';
        pointData3D(nPointProp).data = printProp.potential;
      end
      %
      if isPoromechanics(obj.model)
        obj.VTK.writeVTKFile(printProp.time, pointData3D, cellData3D, [], []);
      elseif isSinglePhaseFlow(obj.model)
        obj.VTK.writeVTKFile(printProp.time, pointData3D, [], [], []);
      end
    end
  end
  
  methods (Static = true)
    % Read the next line and check for eof
    function [flEof,line] = readLine(fid)
      flEof = feof(fid);   % end-of-file flag
      if flEof == 1
        line = '';
      else
        line = strtrim(fgetl(fid));
      end  
    end
  end
end