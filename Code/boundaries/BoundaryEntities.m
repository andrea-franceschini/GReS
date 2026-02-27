classdef BoundaryEntities < handle
  % BOUNDARY CONDITIONS ENTITY class
  % Constrained entities IDs and values of the BCs are defined by the user.

  properties (Access = public)
    % Boundary condition identifier (mainly for error messages)
    name
    % Total number of constrained entities
    totEnts = 0
    % Number of constrained entities for each degree of freedom
    nEntities
    % Indices of constrained entities
    entities
    % Number of input times
    nTimes
    % Set of input times
    times
    % struct with data to set the boundary condition in times
    bcData = struct('time',[],'value',[])
    % Values of currently-stored boundary conditions
    availVals
    % Time id of currently-stored boundary conditions
    availSteps
    % the type of the target entity where the BC is applied
    targetEntity
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
    function obj = BoundaryEntities(name, targetEnt)
      % Calling the function to set object properties
      obj.name = name;
      obj.targetEntity = targetEnt;
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

    function addBCEvent(obj,varargin)

      default = struct('time',[],...
        'value',[]);

      params = readInput(default,varargin{:});

      tVal = params.time;

      if ~isnumeric(tVal)
        time = -1;
      else
        time = tVal;
      end

      obj.times(end+1) = time;
      obj.bcData(end+1) = struct('time', time, 'value', params.value);

      % reorder time in ascending order
      [obj.times,s] = sort(obj.times,"ascend");
      obj.bcData = obj.bcData(s);

      % check for repeated bc times
      if length(unique(obj.times))~=length(obj.times)
        error("Multiple BC events with same time are" + ...
          "not allowed")
      end

    end


  end

  methods (Access = public)

    function setEntities(obj,type,list,comp,mesh)

      if obj.totEnts > 0
        gresLog().warning(2,['Entities for boundary condition %s ' ...
          'already defined. GReS will overwrite them.'])
      end

      [obj.nEntities, obj.entities, obj.entityPos] = readEntitySet(type,list,comp,mesh);

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

    function [nEnts, ents, entsPosition] = readEntitySet(type, ents, components, mesh)
      % read entity set and return the number of entities for each
      % component, the list of entities and their reference location

      switch type
        % input file for list of entities
        case "bclistfile"
          [nEnts, ents] = readListFile(ents);
          entsPosition = getLocation(ents,mesh,obj.obj.targetEntity);
          return
        case {'surfacetags','surfacetag'}
          switch obj.targetEntity
            case "node"
              entsID = unique(mesh.surfaces(ismember(mesh.surfaceTag,surfTags),:));
            case "surface"
              entsID = find(ismember(mesh.surfaceTag,surfTags));
            otherwise
              error("Error for BC %s: XML field surfaceTags is not valid for BC of type %s", obj.name, obj.targetEntity)
          end
        case "bclist"
          entsID = getXMLData(input,[],"bcList");
        case "box"
          boxSize = getXMLData(input,[],"box");
          Lx = boxSize(1:2);
          Ly = boxSize(3:4);
          Lz = boxSize(5:6);
          switch obj.targetEntity
            case "node"
              c = mesh.coordinates;
            case {"cell","volumeforce"}
              c = mesh.cellCentroid;
            otherwise
              error("Error for BC %s: XML field box is not valid for BC of type %s", obj.name, obj.targetEntity)
          end

          entsID = all([ c(:,1) > Lx(1), c(:,1) < Lx(2),...
            c(:,2) > Ly(1), c(:,2) < Ly(2),...
            c(:,3) > Lz(1), c(:,3) < Lz(2)],2);

          entsID = find(entsID);

        otherwise
          error("Unrecognized field 'entityListType' for Boundary condition '%s'.\n" + ...
            "Valid fields are: 'bclist','bclistfile','surfacetags','box'")
      end


      if isempty(entsID)
        error("Error for BC %s: Empty or invalid list of entity in Boundary condition input.",obj.name)
      end

      % expand entity list to components
      compID = true;

      if ~isempty(components)

        if isnumeric(components)
          dir = ["x","y","z"];
          components = dir(components);
        end

        compID =  ismember(["x","y","z"],components);

      end

      if ~any(compID)
        error("Error for BC %s: Check syntax of 'components' field.",obj.name)
      end

      nEnts = numel(entsID).*compID;
      ents = repmat(entsID,sum(compID),1);
      ents = reshape(ents,[],1);

      entsPosition = getLocation(ents,mesh,obj.targetEntity);

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

      switch lower(entType)
        case "node"
          pos = mesh.coordinates(ents,:);
        case "surface"
          pos = mesh.surfaceCentroid(ents,:);
        case {"cell","volumeforce"}
          pos = mesh.cellCentroid(ents,:);
        otherwise
          error("Unrecognized 'targetEntity' for boundary condition '%s':\n" + ...
            "Accepted fields are: 'node','surface','cell','volumeforce'",obj.name)
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
