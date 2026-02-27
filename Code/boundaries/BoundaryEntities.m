classdef BoundaryEntities < handle
  % BOUNDARY CONDITIONS ENTITY class
  % Constrained entities IDs and values of the BCs are defined by the user.

  properties (Access = public)
    % Boundary condition identifier (mainly for error messages)
    name
    % Total number of constrained entities
    totEnts
    % Number of constrained entities for each degree of freedom
    nEntities
    % Indices of constrained entities
    entities
    % Number of input times
    nTimes
    % Set of input times
    times
    % struct with data to set the boundary condition in times
    bcData
    % Values of currently-stored boundary conditions
    availVals
    % Time id of currently-stored boundary conditions
    availSteps
    % the type of the target entity where the BC is applied
    entityType
    % the position in space of the target entities
    entityPos
    % logical index to easily deactivate a bc entity
    isActiveEntity
  end

  properties (Access = private)
    dof
    readSetFlag = true
  end

  methods (Access = public)
    % Class constructor method
    function obj = BoundaryEntities(name, times, bcData, entityType)
      % Calling the function to set object properties
      obj.name = name;
      obj.times = times;
      obj.nTimes = sum(times >= 0.0);
      obj.bcData = bcData;
      obj.entityType = entityType;
    end


    function vals = getValues(obj, t)
      if (obj.nTimes==0)
        % value must a function handle depending on time
        tScale = str2func(['@(t)', char(obj.bcData.time)]);
        if obj.readSetFlag
          obj.availVals(:,1) = readDataSet(obj);
          obj.readSetFlag = false;
        end
        vals = tScale(t)*obj.availVals(:,1); % scale by the time function
      elseif (obj.nTimes == 1)
        if obj.readSetFlag
          obj.availVals(:,1) = readDataSet(obj);
          obj.readSetFlag = false;
        end
        vals = obj.availVals(:,1);
      else
        [i1, i2] = bin_search(obj, t);
        if i1 == i2 % edge case
          if all([i1 i2] == 1)
            obj.availSteps(1) = i1;
            obj.availVals(:,1) = readDataSet(obj,1);
            vals = obj.availVals(:,1);
          elseif all([i1 i2] == obj.nTimes)
            obj.availVals(:,2) = readDataSet(obj,obj.nTimes);
            vals = obj.availVals(:,2);
          end
        else
          p1 = find(obj.availSteps == i1);
          p2 = find(obj.availSteps == i2);
          if (isempty(p1) && isempty(p2))
            p1 = 1;
            obj.availVals(:,1) = readDataSet(obj,i1);
            obj.availSteps(1) = i1;
            p2 = 2;
            obj.availVals(:,2) = readDataSet(obj,i2);
            obj.availSteps(2) = i2;
          elseif (~isempty(p1) && isempty(p2))
            p2 = 3 - p1;
            obj.availVals(:,p2) = readDataSet(obj,i2);
            obj.availSteps(p2) = i2;
          elseif (~isempty(p2) && isempty(p1))
            p1 = 3 - p2;
            obj.availVals(:,p1) = readDataSet(obj,i1);
            obj.availSteps(p1) = i1;
          end
          fac = (t - obj.times(i1)) / (obj.times(i2) - obj.times(i1));
          vals = fac*(obj.availVals(:,p2) - obj.availVals(:,p1)) + obj.availVals(:,p1);
        end

      end

      %vals = vals(obj.isActiveEntity);
    end
  end

  methods (Access = public)

    function setBC(obj, inputStruct, mesh)

      [obj.nEntities, obj.entities, obj.entityPos] = readEntitySet(inputStruct,mesh,obj.entityType,obj.name);

      obj.totEnts = sum(obj.nEntities);

      obj.isActiveEntity = true(obj.totEnts,1);

      if (obj.totEnts == 0)
        error('No boundary conditions are prescribed for %s BC', obj.name);
      end

      obj.availVals = zeros(obj.totEnts,2);
      obj.availSteps = zeros(2,1);
    end

  end

  methods (Access=private)

    function [nEnts, ents, entsPosition] = readEntitySet(type, ents, mesh, entityType, bcName)
      % read entity set (can be a scalar tag or a path to a file)

      if isfield(input,"bcListFile")
        % input list is in a file
        assert(isscalar(fieldnames(input)), ...
          "Error in BC %s: A path to the entityList should be the unique input of the " + ...
          "BCentities field. ", bcName);
        entPath = getXMLData(input,[],"bcListFile");
        [nEnts, ents] = readListFile(entPath);

        entsPosition = getLocation(ents,mesh,entityType);

        return

      else
        entsID = [];
        if isfield(input,"surfaceTags")
          % input is a surface tag with component specification
          surfTags = getXMLData(input,[],"surfaceTags");
          switch entityType
            case "NodeBC"
              entsID = unique(mesh.surfaces(ismember(mesh.surfaceTag,surfTags),:));
            case "SurfBC"
              entsID = find(ismember(mesh.surfaceTag,surfTags));
            otherwise
              error("Error for BC %s: XML field surfaceTags is not valid for BC of type %s", bcName, entityType)
          end
        elseif isfield(input,"bcList")
          % direct entity assignment in the xml file
          entsID = getXMLData(input,[],"bcList");
        elseif isfield(input,"box")
          boxSize = getXMLData(input,[],"box");
          Lx = boxSize(1:2);
          Ly = boxSize(3:4);
          Lz = boxSize(5:6);
          switch entityType
            case "NodeBC"
              c = mesh.coordinates;
            case {"ElementBC","VolumeForce"}
              c = mesh.cellCentroid;
            otherwise
              error("Error for BC %s: XML field box is not valid for BC of type %s", bcName, entityType)
          end

          entsID = all([ c(:,1) > Lx(1), c(:,1) < Lx(2),...
            c(:,2) > Ly(1), c(:,2) < Ly(2),...
            c(:,3) > Lz(1), c(:,3) < Lz(2)],2);

          entsID = find(entsID);
        end

        if isempty(entsID)
          error("Error for BC %s: Empty or invalid list of entity in Boundary condition input.",bcName)
        end

        compID = true;

        if isfield(input,"components")
          dir = getXMLData(input,[],"components");
          compID =  ismember(["x","y","z"],dir);
          if ~any(compID)
            error("Error for BC %s: Check syntax of 'components' field.",bcName)
          end
        end

        nEnts = numel(entsID).*compID;
        ents = repmat(entsID,sum(compID),1);
        ents = reshape(ents,[],1);

        entsPosition = getLocation(ents,mesh,entityType);

      end
    end

    function [nEnts, ents] = readListFile(fileName)

      if (~exist(fileName, 'file'))
        error('File %s does not seem to exist. Please, check the provided file.', fileName);
      end
      header = false;
      fid = fopen(fileName, 'r');
      while (~feof(fid) && ~header)
        line = fgetl(fid);
        word = sscanf(line, '%s');
        if (~strcmp(word(1), '%'))
          % If this is not a commented line (not starting with %)
          nEnts = sscanf(line, '%i');
          header = true;
        end
      end
      if (~header)
        error('Missing header in readSet.');
      end
      nValMax = sum(nEnts);
      ents = zeros(nValMax,1);
      id = 1;
      while ~feof(fid)
        line = fgetl(fid);
        word = sscanf(line, '%s');
        if (~strcmp(word(1), '%'))
          % If this is not a commented line (not starting with %)
          num = sscanf(line, '%i');
          nNum = length(num);
          ents(id:id+nNum-1) = num;
          id = id + nNum;
        end
      end
      fclose(fid);

    end


    function pos = getLocation(ents,mesh,entType)

      switch entType
        case "NodeBC"
          pos = mesh.coordinates(ents,:);
        case "SurfBC"
          pos = mesh.surfaceCentroid(ents,:);
        case {"ElementBC","VolumeForce"}
          pos = mesh.cellCentroid(ents,:);
      end
    end

    function [i1, i2] = bin_search(obj, t)
      % Return the interval [i1, i2] such that obj.times(i1) <= t <= obj.times(i2)

      if t <= obj.times(1)
        i1 = 1;
        i2 = 1;
        return;
      elseif t >= obj.times(end)
        i1 = obj.nTimes;
        i2 = obj.nTimes;
        return;
      end

      i1 = 1;
      i2 = obj.nTimes;

      while (i2 - i1 > 1)
        pos = floor((i1 + i2)/2);
        if obj.times(pos) > t
          i2 = pos;
        else
          i1 = pos;
        end
      end
    end

  end
end
