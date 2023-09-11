classdef OutState < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = public)
    timeList
    m
  end
  
  properties (Access = private)
    model
    printPath
    mesh
    timeID = 1
    VTK
    flOutData = true
    material
%     flPrint = true
  end
  
  methods (Access = public)
    function obj = OutState(symMod,mat,grid,fileName,varargin)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      obj.setOutState(symMod,mat,grid,fileName,varargin)
    end
    
    function finalize(obj)
      obj.VTK.finalize();
    end
    
    function printState(obj,stateOld,stateNew)
      if nargin == 2
        printProp.time = stateOld.t;
%         if isPoromechanics(obj.model) && isSinglePhaseFlow(obj.model)
%           [avStressOld,avStrainOld] = finalizeStatePoro(stateOld);
%           [fluidPotOld] = finalizeStateFlow(stateOld);
%           printProp = struct('time',stateOld.t,'displ',stateOld.displ, ...
%             'stress',avStressOld,'strain',avStrainOld, ...
%             'pressure',stateOld.pressure,'potential',fluidPotOld);
        if obj.flOutData
          obj.m.expTime(obj.timeID,1) = printProp.time;
        end
        if isPoromechanics(obj.model)
          [avStressOld,avStrainOld] = finalizeStatePoro(stateOld);
%           printProp = struct('time',stateOld.t,'displ',stateOld.displ, ...
%             'stress',avStressOld,'strain',avStrainOld);
          printProp.dispConv = stateOld.dispConv;
          printProp.stress = avStressOld;
          printProp.strain = avStrainOld;
          if obj.flOutData
            obj.m.expDispl(:,obj.timeID) = printProp.dispConv;
          end
        end
        if isFlow(obj.model)
          [fluidPotOld] = finalizeStateFlow(stateOld);
%           printProp = struct('time',stateOld.t,'pressure',stateOld.pressure, ...
%           'potential',fluidPotOld);
          printProp.pressure = stateOld.pressure;
          printProp.potential = fluidPotOld;
          if isVariabSatFlow(obj.model)
            printProp.Sw = stateOld.watSat;
          end
          if obj.flOutData
            obj.m.expPress(:,obj.timeID) = printProp.pressure;
            if isVariabSatFlow(obj.model)
              obj.m.expSw(:,obj.timeID) = printProp.Sw;
            end
          end
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
            printProp.time = obj.timeList(obj.timeID);
%             if isPoromechanics(obj.model) && isSinglePhaseFlow(obj.model)
%               [avStressOld,avStrainOld] = finalizeStatePoro(stateOld);
%               [avStressNew,avStrainNew] = finalizeStatePoro(stateNew);
%               [fluidPotOld] = finalizeStateFlow(stateOld);
%               [fluidPotNew] = finalizeStateFlow(stateNew);
%               printProp = struct('time',obj.timeList(obj.timeID), ...
%                 'displ',stateNew.displ*fac+stateOld.displ*(1-fac), ...
%                 'stress',avStressNew*fac+avStressOld*(1-fac), ...
%                 'strain',avStrainNew*fac+avStrainOld*(1-fac), ...
%                 'pressure',stateNew.pressure*fac+stateOld.pressure*(1-fac), ...
%                 'potential',fluidPotNew*fac+fluidPotOld*(1-fac));
            if obj.flOutData
              obj.m.expTime(obj.timeID+1,1) = printProp.time;
            end
            if isPoromechanics(obj.model)
              [avStressOld,avStrainOld] = finalizeStatePoro(stateOld);
              [avStressNew,avStrainNew] = finalizeStatePoro(stateNew);
%               printProp = struct('time',obj.timeList(obj.timeID), ...
%                 'displ',stateNew.displ*fac+stateOld.displ*(1-fac), ...
%                 'stress',avStressNew*fac+avStressOld*(1-fac), ...
%                 'strain',avStrainNew*fac+avStrainOld*(1-fac));
              printProp.dispConv = stateNew.dispConv*fac+stateOld.dispConv*(1-fac);
              printProp.stress = avStressNew*fac+avStressOld*(1-fac);
              printProp.strain = avStrainNew*fac+avStrainOld*(1-fac);
              if obj.flOutData
                obj.m.expDispl(:,obj.timeID+1) = printProp.dispConv;
              end
            end
            if isFlow(obj.model)
              [fluidPotOld] = finalizeStateFlow(stateOld);
              [fluidPotNew] = finalizeStateFlow(stateNew);
%               printProp = struct('time',obj.timeList(obj.timeID), ...
%                 'pressure',stateNew.pressure*fac+stateOld.pressure*(1-fac), ...
%                 'potential',fluidPotNew*fac+fluidPotOld*(1-fac));
              printProp.pressure = stateNew.pressure*fac+stateOld.pressure*(1-fac);
              printProp.potential = fluidPotNew*fac+fluidPotOld*(1-fac);
              if obj.flOutData
                obj.m.expPress(:,obj.timeID+1) = printProp.pressure;
              end
              if isVariabSatFlow(obj.model)
%                 printProp.Sw = stateNew.watSat*fac+stateOld.watSat*(1-fac);
                printProp.Sw = obj.material.computeSwAnddSw(obj.mesh,printProp.pressure);
                if obj.flOutData
                  obj.m.expSw(:,obj.timeID+1) = printProp.Sw;
                end
              end
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
    function setOutState(obj,symMod,mat,grid,fileName,data)
      if ~isempty(data)
        obj.GaussPts = data{1};
      end
      obj.model = symMod;
      obj.material = mat;
      obj.mesh = grid.topology;
      %
      fid = fopen(fileName,'r');
      [flEof,line] = OutState.readLine(fid);
      %
      block = '';
      while ~strcmp(line,'End')
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
        obj.timeList = [];
      end
      fclose(fid);
      %
      obj.VTK = VTKOutput(obj.mesh);
      %
      % TO DO: Here we are printing everything; consider adding the
      % possibility of printing only the results of a portion of the domain
      if obj.flOutData
        if isfile('expData.mat')
          delete 'expData.mat'
        end
        obj.m = matfile('expData.mat','Writable',true);
        l = length(obj.timeList) + 1;
        obj.m.expTime = zeros(l,1);
        if isFlow(obj.model)
          if isFEMBased(obj.model,'Flow')
            obj.m.expPress = zeros(obj.mesh.nNodes,l);
          elseif isFVTPFABased(obj.model,'Flow')
            obj.m.expPress = zeros(obj.mesh.nCells,l);
            if isVariabSatFlow(obj.model)
              obj.m.expSw = zeros(obj.mesh.nCells,l);
            end
          end
        end
        if isPoromechanics(obj.model)
          obj.m.expDispl = zeros(obj.mesh.nDim*obj.mesh.nNodes,l);
          % Maybe consider adding other output properties
        end
      end
    end
    %
    
    function buildPrintStruct(obj,printProp)
      nPointProp = 0;
      nCellProp = 0;
%       if isCoupFlowPoro(obj.model)
%         % Poromechanics
%         nPointProp = 3;
%         nCellProp = 12;
%         % Flow part
%         if isFEMBased(obj.model)
%           nPointProp = nPointProp + 2;
%         elseif isFVTPFABased(obj.model)
%           nCellProp = nCellProp + 2;
%         end
      if isPoromechanics(obj.model)
        nPointProp = nPointProp + 3;
        nCellProp = nCellProp + 12;
      end
      if isFlow(obj.model)
        if isFEMBased(obj.model,'Flow')
          nPointProp = nPointProp + 2;   %2
        elseif isFVTPFABased(obj.model,'Flow')
          nCellProp = nCellProp + 2;
        end
        if isVariabSatFlow(obj.model)
          nCellProp = nCellProp + 1;
        end
      end
      %
      if nPointProp > 0
        pointData3D = repmat(struct('name', 1, 'data', 1), nPointProp, 1);
      else
        pointData3D = [];
      end
      if nCellProp > 0
        cellData3D = repmat(struct('name', 1, 'data', 1), nCellProp, 1);
      else
        cellData3D = [];
      end
      %
      if isPoromechanics(obj.model)
        %
        % Displacement
        pointData3D(1).name = 'ux';
        pointData3D(1).data = printProp.dispConv(1:3:length(printProp.dispConv));
        pointData3D(2).name = 'uy';
        pointData3D(2).data = printProp.dispConv(2:3:length(printProp.dispConv));
        pointData3D(3).name = 'uz';
        pointData3D(3).data = printProp.dispConv(3:3:length(printProp.dispConv));
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
      if isFlow(obj.model)
        %
        % Pressure and potential
        if isFEMBased(obj.model,'Flow')
          pointData3D(nPointProp-1).name = 'press';   %-1
          pointData3D(nPointProp-1).data = printProp.pressure;
          pointData3D(nPointProp).name = 'potential';
          pointData3D(nPointProp).data = printProp.potential;
%           pointData3D(nPointProp).name = 'nodeID';
%           pointData3D(nPointProp).data = (1:length(printProp.pressure))';
        elseif isFVTPFABased(obj.model,'Flow')
          cellData3D(nCellProp-1).name = 'press';
          cellData3D(nCellProp-1).data = printProp.pressure;
          cellData3D(nCellProp).name = 'potential';
          cellData3D(nCellProp).data = printProp.potential;
          if isVariabSatFlow(obj.model)
            cellData3D(nCellProp-2).name = 'Sw';
            cellData3D(nCellProp-2).data = printProp.Sw;
          end
        end
      end
      %
%       if isPoromechanics(obj.model)
        obj.VTK.writeVTKFile(printProp.time, pointData3D, cellData3D, [], []);
%       elseif isSinglePhaseFlow(obj.model)
%         obj.VTK.writeVTKFile(printProp.time, pointData3D, [], [], []);
%       end
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