classdef OutState_new < handle
   % Class for printing results to VTK
   % Input: OutState(model,mesh,fNameList)
   % Optional parameters:
   % 'writeMat',true/false: print result to mat-file
   % 'outFolder','folderName': specify name of output folder



   properties (Access = public)
      modTime = false %flag for time step size matching timeList
      timeList
      results
      model
      writeSolution
   end

   properties (Access = private)
      timeID = 1
      VTK
   end

   methods (Access = public)
      function obj = OutState_new(model,mesh,fileName,options)
         arguments
            model (1,1) ModelType
            mesh (1,1) Mesh
            fileName (1,1) string
            options.flagMatFile logical = true
            options.folderName string = "vtkOutput"
         end
         % deal variable input
         obj.writeSolution = options.flagMatFile;
         foldName = options.folderName;
         obj.setOutState(model,mesh,fileName,foldName)
      end

      function finalize(obj)
         obj.VTK.finalize();
      end

      function printState_new(obj,solv,stateOld,stateNew)
         % print solution of the model according to the print time in the
         % list
         % Initialize input structure for VTK output
         cellData3D = [];
         pointData3D = [];
         if nargin == 3
            time = stateOld.t;
            % print result to mat-file (to be updated in future version of
            % the code)
            if obj.writeSolution
               obj.results.expTime(obj.timeID,1) = time;
               if isPoromechanics(obj.model)
                  obj.results.expDispl(:,obj.timeID) = stateOld.dispConv;
               end
               if isFlow(obj.model)
                  obj.results.expPress(:,obj.timeID) = stateOld.pressure;
               end
               if isVariabSatFlow(obj.model)
                  obj.results.expSw(:,obj.timeID) = stateOld.Sw;
               end
            end
            for fld = solv.fields
               [cellData,pointData] = printState(solv.getSolver(fld),stateOld);
               cellData3D = [cellData3D; cellData];
               pointData3D = [pointData3D; pointData];
            end
            obj.VTK.writeVTKFile(time, pointData3D, cellData3D, [], []);
            % update the print structure
         elseif nargin == 4
            if obj.timeID <= length(obj.timeList)
               while (obj.timeList(obj.timeID) <= stateNew.t)
                  assert(obj.timeList(obj.timeID) > stateOld.t, ...
                     'Print time %f out of range (%f - %f)',obj.timeList(obj.timeID), ...
                     stateOld.t,stateNew.t);
                  assert(stateNew.t - stateOld.t > eps('double'),'Dt too small for printing purposes');
                  %
                  time = obj.timeList(obj.timeID);
                  if obj.writeSolution
                     % print solution to mat-file
                     fac = (time - stateOld.t)/(stateNew.t - stateOld.t);
                     obj.results.expTime(obj.timeID+1,1) = time;
                     if isPoromechanics(obj.model)
                        obj.results.expDispl(:,obj.timeID+1) = stateNew.dispConv*fac+stateOld.dispConv*(1-fac);
                     end
                     if isFlow(obj.model)
                        obj.results.expPress(:,obj.timeID+1) = stateNew.pressure*fac+stateOld.pressure*(1-fac);
                     end
                     if isVariabSatFlow(obj.model)
                        obj.results.expSw(:,obj.timeID) = stateNew.watSat*fac+stateOld.watSat*(1-fac);
                     end
                  end
                  % Write output structure looping trough available models
                  for fld = solv.fields
                     [cellData,pointData] = printState(solv.getSolver(fld),...
                        stateOld, stateNew, time);
                     cellData3D = [cellData3D; cellData];
                     pointData3D = [pointData3D; pointData];
                  end
                  obj.VTK.writeVTKFile(time, pointData3D, cellData3D, [], []);
                  obj.timeID = obj.timeID + 1;
                  if obj.timeID > length(obj.timeList)
                     break
                  end
               end
            end
         end
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
                if obj.writeSolution
                    obj.results.expTime(obj.timeID,1) = printProp.time;
                end
                if isPoromechanics(obj.model)
                    [avStressOld,avStrainOld] = finalizeStatePoro(stateOld);
                    %           printProp = struct('time',stateOld.t,'displ',stateOld.displ, ...
                    %             'stress',avStressOld,'strain',avStrainOld);
                    printProp.dispConv = stateOld.dispConv;
                    printProp.stress = avStressOld;
                    printProp.strain = avStrainOld;
                    if obj.writeSolution
                        obj.results.expDispl(:,obj.timeID) = printProp.dispConv;
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
                    if obj.writeSolution
                        obj.results.expPress(:,obj.timeID) = printProp.pressure;
                        if isVariabSatFlow(obj.model)
                            obj.results.expSw(:,obj.timeID) = printProp.Sw;
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
                        if obj.writeSolution
                            obj.results.expTime(obj.timeID+1,1) = printProp.time;
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
                            if obj.writeSolution
                                obj.results.expDispl(:,obj.timeID+1) = printProp.dispConv;
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
                            if obj.writeSolution
                                obj.results.expPress(:,obj.timeID+1) = printProp.pressure;
                            end
                            if isVariabSatFlow(obj.model)
                                %                 printProp.Sw = stateNew.watSat*fac+stateOld.watSat*(1-fac);
                                printProp.Sw = obj.material.computeSwAnddSw(obj.mesh,printProp.pressure);
                                if obj.writeSolution
                                    obj.results.expSw(:,obj.timeID+1) = printProp.Sw;
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
        function setOutState(obj,model,mesh,fileName,foldName)
           obj.model = model;
          %
          readTimeList(obj,fileName);
          %
          obj.VTK = VTKOutput(mesh,foldName);
          % Write solution to matfile. This feature will be extended in a
          % future version of the code
          if obj.writeSolution
             if isfile('expData.mat')
                delete 'expData.mat'
             end
             obj.results = matfile('expData.mat','Writable',true);
             setMatFile(obj,mesh);
          end
       end
       %


       function buildPrintStruct(obj,printProp)
         
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

       function readTimeList(obj,fileName)
          fid = fopen(fileName,'r');
          [flEof,line] = OutState.readLine(fid);
          if isempty(sscanf(line,'%f'))
             line = strtok(line);
             if strcmpi(line,'on')
                obj.modTime = true;
             elseif strcmpi(line,'off')
                obj.modTime = false;
             end
             if ~strcmpi(line,'end')
                [flEof,line] = OutState.readLine(fid);
             end
          end

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
       end

       function setMatFile(obj,msh)
          l = length(obj.timeList) + 1;
          obj.results.expTime = zeros(l,1);
          if isFlow(obj.model)
             if isFEMBased(obj.model,'Flow')
                obj.results.expPress = zeros(msh.nNodes,l);
             elseif isFVTPFABased(obj.model,'Flow')
                obj.results.expPress = zeros(msh.nCells,l);
                if isVariabSatFlow(obj.model)
                   obj.results.expSw = zeros(msh.nCells,l);
                end
             end
          end
          if isPoromechanics(obj.model)
             obj.results.expDispl = zeros(msh.nDim*msh.nNodes,l);
             % Maybe consider adding other output properties
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