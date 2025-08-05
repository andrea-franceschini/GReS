classdef OutState < handle
   % Class for printing results to VTK
   % Input: OutState(model,mesh,fNameList)
   % Optional parameters:
   % 'flagMatFile',true/false: print result to mat-file
   % 'folderName','folderName': specify name of output folder

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
      writeVtk
   end

   methods (Access = public)
      function obj = OutState(model,mesh,fileName,options)
         arguments
            model (1,1) ModelType
            mesh (1,1) Mesh
            fileName (1,1) string
            options.flagMatFile logical = true
            options.writeVtk logical = true
            options.folderName string = "vtkOutput"
         end
         % deal variable input
         obj.writeVtk = options.writeVtk;
         obj.writeSolution = options.flagMatFile;
         foldName = options.folderName;
         obj.setOutState(model,mesh,fileName,foldName);
      end

      function finalize(obj)
          if obj.writeVtk
            obj.VTK.finalize();
          end
      end

      function printState(obj,solv,bound,stateOld,stateNew)
         % print solution of the model according to the print time in the
         % list
         % Initialize input structure for VTK output
         cellData3D = [];
         pointData3D = [];
         if nargin == 4
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
                  obj.results.expSat(:,obj.timeID) = stateOld.saturation;
               end
            end
            fields = [bound.field];
            for fld = solv.fields
               idx = find(fields == fld);
               [cellData,pointData] = printState(solv.getSolver(fld),bound(idx),stateOld);
               cellData3D = [cellData3D; cellData];
               pointData3D = [pointData3D; pointData];
            end
            if obj.writeVtk
               obj.VTK.writeVTKFile(time, pointData3D, cellData3D, [], []);
            end
            
            % update the print structure
         elseif nargin == 5
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
                     % obj.results.expTime(obj.timeID+1,1) = time;
                     obj.results(obj.timeID+1).expTime = time;
                     if isPoromechanics(obj.model)
                        % obj.results.expDispl(:,obj.timeID+1) = stateNew.dispConv*fac+stateOld.dispConv*(1-fac);
                        obj.results(obj.timeID+1).expDispl = stateNew.dispConv*fac+stateOld.dispConv*(1-fac);
                     end
                     if isFlow(obj.model)
                        % obj.results.expPress(:,obj.timeID+1) = stateNew.pressure*fac+stateOld.pressure*(1-fac);
                        obj.results(obj.timeID+1).expPress = stateNew.pressure*fac+stateOld.pressure*(1-fac);
                     end
                     if isVariabSatFlow(obj.model)
                        % obj.results.expSat(:,obj.timeID+1) = stateNew.saturation*fac+stateOld.saturation*(1-fac);
                        obj.results(obj.timeID+1).expSat = stateNew.saturation*fac+stateOld.saturation*(1-fac);
                     end
                  end
                  % Write output structure looping trough available models
                  for fld = solv.fields                     
                     [cellData,pointData] = printState(solv.getSolver(fld),...
                        bound, stateOld, stateNew, time);
                     % merge new fields
                     cellData3D = OutState.mergeOutFields(cellData3D,cellData);
                     pointData3D = OutState.mergeOutFields(pointData3D,pointData);
                  end
                  if obj.writeVtk
                    obj.VTK.writeVTKFile(time, pointData3D, cellData3D, [], []);
                  end
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
         obj.timeList = OutState.readTime(fileName);
         %
         if obj.writeVtk
            obj.VTK = VTKOutput(mesh,foldName);
         end
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

      function tList = readTimeList(obj,fileName)
         fid = fopen(fileName,'r');
         [flEof,line] = OutState.readLine(fid);
         if isempty(sscanf(line,'%f'))
            line = strtok(line);
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
            tList = str2double(blockSplt);
            if any(isnan(obj.timeList))
               error('There are invalid entries in the list of output times')
            end
         else
            tList = [];
         end
         fclose(fid);
      end

      function setMatFile(obj,msh)
         l = length(obj.timeList) + 1;
         obj.results = repmat(struct('expTime', 0),l,1);
         for i=1:l
            if isFlow(obj.model)
               obj.results(i).expPress=[];
               if isVariabSatFlow(obj.model)
                  obj.results(i).expSat=[];
               end
            end
            if isPoromechanics(obj.model)
               obj.results(i).expDispl = [];
            end
         end

         % obj.results.expTime = zeros(l,1);
         % if isFlow(obj.model)
         %    if isFEMBased(obj.model,'Flow')
         %       obj.results.expPress = zeros(msh.nNodes,l);
         %    elseif isFVTPFABased(obj.model,'Flow')
         %       obj.results.expPress = zeros(msh.nCells,l);
         %       if isVariabSatFlow(obj.model)
         %          obj.results.expSat = zeros(msh.nCells,l);
         %       end
         %    end
         % end
         % if isPoromechanics(obj.model)
         %    obj.results.expDispl = zeros(msh.nDim*msh.nNodes,l);
         %    % Consider adding other output properties
         % end
      end
   end

   methods (Static = true)

       function mergeStruct = mergeOutFields(strA,strB)
           % Merge output variable coming from a field to the global output
           % structure
           % Concatenate the two structure arrays
           if isempty(strA)
               mergeStruct = strB;
               return
           end
           mergeStruct = [strA; strB];
           names = {mergeStruct.name};
           [~, uniqueIdx] = unique(names, 'stable');
           % Create the merged structure using the unique indices
           mergeStruct = mergeStruct(uniqueIdx);
       end

       function [flEof,line] = readLine(fid)
           % Read the next line and check for eof
           flEof = feof(fid);   % end-of-file flag
           if flEof == 1
               line = '';
           else
               line = strtrim(fgetl(fid));
           end
       end

       function tList = readTime(fileName)
           fid = fopen(fileName,'r');
           [flEof,line] = OutState.readLine(fid);
           if isempty(sscanf(line,'%f'))
               line = strtok(line);
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
               tList = str2double(blockSplt);
               if any(isnan(tList))
                   error('There are invalid entries in the list of output times')
               end
           else
               tList = [];
           end
           fclose(fid);
       end
   end
end