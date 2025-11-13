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
    timeID = 1
    VTK
    writeVtk
  end

  properties (Access = private)
    % timeID = 1
    % VTK
    % writeVtk
  end

  methods (Access = public)

    function obj = OutState(model, mesh, varargin)

      % Default parameters
      flagMatFile = false;
      writeVtk = true;
      folderName = "vtkOutput";
      tList = [];

      % ------------------------------------------------------------
      % Detect input mode
      % ------------------------------------------------------------
      if isempty(varargin)
        % no output
        obj.writeVtk = false;
        obj.writeSolution = false;
        obj.model = model;
        return
      end

      firstArg = varargin{1};

      if isscalar(varargin) && (ischar(firstArg) || isstring(firstArg) || isstruct(firstArg))
        % ---------------- XML configuration mode ----------------
        xml = firstArg;
        if ~isstruct(xml)
          xml = readstruct(xml, AttributeSuffix="");
        end

        flagMatFile = logical(getXMLData(xml, 0, "saveHistory"));
        folderName  = getXMLData(xml, "vtkOutput", "outputFile");
        tList       = getXMLData(xml, [], "printTimes");

        % If no outputFile provided, disable VTK
        writeVtk = isfield(xml, "outputFile");

      else
        % ---------------- Key-value pair mode ----------------
        if mod(numel(varargin), 2) ~= 0
          error('OutState:InvalidInput', ...
            'Key-value pair inputs must come in pairs.');
        end

        opt = struct(varargin{:});
        if isfield(opt, 'flagMatFile'), flagMatFile = logical(opt.flagMatFile); end
        if isfield(opt, 'writeVtk'),    writeVtk    = logical(opt.writeVtk);    end
        if isfield(opt, 'folderName'),  folderName  = string(opt.folderName);   end
        if isfield(opt, 'timeList'),    tList       = opt.timeList;             end
      end

      % ------------------------------------------------------------
      % Object setup
      % ------------------------------------------------------------
      obj.writeVtk = writeVtk;
      obj.writeSolution = flagMatFile;
      obj.model = model;
      obj.VTK = VTKOutput(mesh, folderName);

      % MAT-file output
      if obj.writeSolution
        if isfile('expData.mat'), delete 'expData.mat'; end
        obj.results = matfile('expData.mat', 'Writable', true);
        setMatFile(obj, mesh);
      end

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

    end


    function finalize(obj)
      if obj.writeVtk
        obj.VTK.finalize();
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

    function setMatFile(obj,msh)
      l = length(obj.timeList) + 1;
      obj.results = repmat(struct('expTime', 0),l,1);
      for i=1:l
        if isFlow(obj.model)
          if isFEMBased(obj.model,'Flow')
            obj.results(i).expPress = [];
          elseif isFVTPFABased(obj.model,'Flow')
            obj.results(i).expPress = [];
            if isVariabSatFlow(obj.model)
              obj.results(i).expSat = [];
            end
          end
        end
        if isPoromechanics(obj.model)
          obj.results(i).expDispl = [];
          % Consider adding other output properties
        end
      end
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