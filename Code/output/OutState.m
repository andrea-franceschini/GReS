classdef OutState < handle & matlab.mixin.Copyable
  % Class for printing results to VTK
  % Input: OutState(model,mesh,fNameList)
  % Optional parameters:
  % 'flagMatFile',true/false: print result to mat-file
  % 'folderName','folderName': specify name of output folder

  properties (Access = public)
    matFile
    % model
    timeID = 1
    vtkFile
  end

  properties (SetAccess = private)
    writeVtk
    timeList
    matFileName
    vtkFileName
    writeSolution
    isFolderReady = false
  end

  methods (Access = public)

    function obj = OutState(varargin)

      if nargin == 0
        return
      end

      % Default parameters
      folderName = "output";
      tList = [];

      % ------------------------------------------------------------
      % Detect input mode
      % ------------------------------------------------------------
      if isempty(varargin)
        % no output
        obj.writeVtk = false;
        obj.writeSolution = false;
        % obj.model = model;
        return
      end

      firstArg = varargin{1};

      if isscalar(varargin) && (ischar(firstArg) || isstring(firstArg) || isstruct(firstArg))
        % ---------------- XML configuration mode ----------------
        input = firstArg;
        if ~isstruct(input)
          input = readstruct(input, AttributeSuffix="");
        end

        if isfield(input,"Output")
          input = input.Output;
        end

        if isfield(input,"fileName")
          assert(isscalar(fieldnames(input)),"FileName, " + ...
            " must be a unique parameter.");
          input = readstruct(input.fileName,AttributeSuffix="");
        end

        obj.writeSolution = logical(getXMLData(input, 0, "saveHistory"));
        obj.vtkFileName  = getXMLData(input, folderName, "outputFile");
        obj.matFileName = getXMLData(input, folderName, "matFileName");

        % If no outputFile name provided, disable VTK
        obj.writeVtk = isfield(input, "outputFile");

        if any([obj.writeSolution,obj.writeVtk])
          tList = getXMLData(input, [], "printTimes");
        else
          return
        end


      else
        % ---------------- Key-value pair mode ----------------
        if mod(numel(varargin), 2) ~= 0
          error('OutState:InvalidInput', ...
            'Key-value pair inputs must come in pairs.');
        end

        opt = struct(varargin{:});
        if isfield(opt, 'flagMatFile'), obj.writeSolution = logical(opt.flagMatFile); end
        if isfield(opt, 'writeVtk'),    obj.writeVtk    = logical(opt.writeVtk);    end
        if isfield(opt, 'folderName'),  obj.vtkFileName  = string(opt.folderName);   end
        if isfield(opt, 'matFileName'),  obj.matFileName  = string(opt.matFileName);   end
        if isfield(opt, 'timeList'),    tList       = opt.timeList;             end
      end

      % ------------------------------------------------------------
      % Object setup
      % ------------------------------------------------------------

      % Time list handling
      if obj.writeVtk
        % vtm file document node
        assert(~isempty(tList), ...
          "Print times have not been specified for output.");
      end

      if isnumeric(tList)
        obj.timeList = reshape(tList, [], 1);
      else
        obj.timeList = obj.readTimeList(tList);
      end

      % MAT-file output
      if obj.writeSolution
        if isfile('expData.mat'), delete 'expData.mat'; end
        nT = length(obj.timeList);
        obj.matFile = repmat(struct('time', 0),nT,1);
      end

    end

    function prepareOutputFolders(obj)

      foldName = sprintf('%s/output_%5.5i',obj.vtkFileName,obj.timeID);

      if (~obj.isFolderReady)
        createVTKFolder(obj);
      end

      % create folder for vtu files
      status = mkdir(foldName);
      if (status ~= 1)
        error('Unable to create folder for VTK output.');
      end
    end


    function writeVTKfile(obj,block,vtuName,mesh,time,pointData3D,cellData3D,pointData2D,cellData2D)

      % block: the xml block of the vtm file in which we write the vtu dataset
      % data struct: struct array with fields 'name' and 'data'

      vtmName = sprintf('output_%5.5i/%s',obj.timeID,vtuName);
      outName = sprintf('%s/%s.vtu',obj.vtkFileName,vtmName);

      % call mex vtk writer
      if ~all(isempty([cellData3D; pointData3D]))
        mxVTKWriter(outName, time, mesh.coordinates, mesh.cells, mesh.cellVTKType, ...
          mesh.cellNumVerts, pointData3D, cellData3D);
      elseif ~all(isempty([cellData2D; pointData2D]))
        mxVTKWriter(outName, time, mesh.coordinates, mesh.surfaces, mesh.surfaceVTKType, ...
          mesh.surfaceNumVerts, pointData2D, cellData2D);
      end

      % write dataset to vtm block
      dataset = obj.vtkFile.createElement('DataSet');
      dataset.setAttribute('name', vtuName);
      dataset.setAttribute('file', sprintf('%s.vtu',vtmName));
      % append dataset to vtm block
      block.appendChild(dataset);
    end

    function writeVTMFile(obj)
      fname = sprintf('%s/output_%5.5i.vtm',obj.vtkFileName,obj.timeID);
      xmlwrite(fname, obj.vtkFile);

    end

    function finalize(obj)

      % write the pvd file
      pvd = com.mathworks.xml.XMLUtils.createDocument('VTKFile');
      toc = pvd.getDocumentElement;

      toc.setAttribute('type', 'Collection');
      toc.setAttribute('version', '1.0');
      blocks = pvd.createElement('Collection');


      for i = 1 : obj.timeID-1
        block = pvd.createElement('DataSet');
        block.setAttribute('timestep', sprintf('%e', obj.timeList(i)));
        [~,fname,~] = fileparts(obj.vtkFileName);
        % standard naming for vtm files
        vtmFileName = sprintf('%s/output_%5.5i.vtm',fname,i);
        block.setAttribute('file', vtmFileName);
        blocks.appendChild(block);
      end

      toc.appendChild(blocks);

      fileName = sprintf('%s.pvd', obj.vtkFileName);
      xmlwrite(fileName, pvd);

      if obj.writeSolution
        output = obj.matFile;
        save(strcat(obj.matFileName,'.mat'),"output")
      end

    end

  end

  methods (Access = private)

   
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

    function createVTKFolder(obj)
      if (isfolder(obj.vtkFileName))
        rmdir(obj.vtkFileName, 's');
      end

      status = mkdir(obj.vtkFileName);
      if (status ~= 1)
        error('Unable to create folder for VTK output.');
      end
      obj.isFolderReady = true;
    end



  end

  methods (Static = true)

    function mergeStruct = mergeOutFields(strA,strB)
      % Merge output variable coming from a field to the global output
      % structure
      % Concatenate the two structure arrays
      strA = reshape(strA,[],1);
      strB = reshape(strB,[],1);

      if isempty(strA)
        mergeStruct = strB;
      elseif isempty(strB)
        mergeStruct = strA;
      else

        try
          mergeStruct = [strA; strB];
        catch
          unknownFields = setdiff(fieldnames(strB),fieldnames(strA));
          error("Ouput structure must be an array of structures with 2 fields:" + ...
            " 'name' and 'data'. \n" + ...
            "Invalid fields: \n %s, ",unknownFields{:})

        end
      end

      names = {mergeStruct.name};

      % remove empty 
      isEmpty = cellfun(@isempty, names);
      mergeStruct(isEmpty) = [];
      names(isEmpty) = [];

      [~, uniqueIdx] = unique(names, 'stable');
      % Create the merged structure using the unique indices
      mergeStruct = mergeStruct(uniqueIdx);
    end

    function outData = printMeshData(mesh,data)
      cellStr = repmat(struct('name', 1, 'data', 1), 1, 1);
      % Displacement
      cellStr(1).name = 'cellTag';
      cellStr(1).data = mesh.cellTag;
      outData = OutState.mergeOutFields(data,cellStr);
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