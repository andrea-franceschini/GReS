classdef OutState < handle
  % Class for printing results to VTK
  % Input: OutState(model,mesh,fNameList)
  % Optional parameters:
  % 'flagMatFile',true/false: print result to mat-file
  % 'folderName','folderName': specify name of output folder

  properties (Access = public)
    modTime = false     %flag for time step size matching timeList
    timeList
    results
    % model
    writeSolution
    timeID = 1
    VTK
    writeVtk
    matFileName
    vtkFileName
  end

  properties (Access = private)
    % timeID = 1
    % VTK
    % writeVtk
  end

  methods (Access = public)

    function obj = OutState(mesh, varargin)

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

      %obj.model = model;
      obj.VTK = VTKOutput(mesh, obj.vtkFileName);

      % Time list handling
      if obj.writeVtk
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
        %obj.results = matfile('expData.mat', 'Writable', true);
        nT = length(obj.timeList);
        obj.results = repmat(struct('time', 0),nT,1);
      end

    end


    function finalize(obj)
      if obj.writeVtk
        obj.VTK.finalize();
      end

      if obj.writeSolution
        output = obj.results;
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

    % move this logic into specific Physics solver
    % function setMatFile(obj,msh)
    %   l = length(obj.timeList) + 1;
    %   obj.results.expTime = zeros(l,1);
    %   if isFlow(obj.model)
    %     if isFEMBased(obj.model,'Flow')
    %       obj.results.expPress = zeros(msh.nNodes,l);
    %     elseif isFVTPFABased(obj.model,'Flow')
    %       obj.results.expPress = zeros(msh.nCells,l);
    %       if isVariabSatFlow(obj.model)
    %         obj.results.expSat = zeros(msh.nCells,l);
    %       end
    %     end
    %   end
    %   if isPoromechanics(obj.model)
    %     obj.results.expDispl = zeros(msh.nDim*msh.nNodes,l);
    %     % Consider adding other output properties
    %   end
    % end

    % function setMatFile(obj,msh)
    %   l = length(obj.timeList) + 1;
    %   obj.results = repmat(struct('expTime', 0),l,1);
    %   for i=1:l
    %     if isFlow(obj.model)
    %       if isFEMBased(obj.model,'Flow')
    %         obj.results(i).expPress = [];
    %       elseif isFVTPFABased(obj.model,'Flow')
    %         obj.results(i).expPress = [];
    %         if isVariabSatFlow(obj.model)
    %           obj.results(i).expSat = [];
    %         end
    %       end
    %     end
    %     if isPoromechanics(obj.model)
    %       obj.results(i).expDispl = [];
    %       % Consider adding other output properties
    %     end
    %   end
    % end

  end

  methods (Static = true)

    function mergeStruct = mergeOutFields(strA,strB)
      % Merge output variable coming from a field to the global output
      % structure
      % Concatenate the two structure arrays

      try
        mergeStruct = [strA; strB];
      catch
        unknownFields = setdiff(fieldnames(strB),fieldnames(strA));
        error("Ouput structure must be an array of structures with 2 fields:" + ...
          " 'name' and 'data'. \n" + ...
          "Invalid fields: \n %s, ",unknownFields{:})

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