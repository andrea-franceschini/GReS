classdef BoundaryEntities < handle
  % BOUNDARY CONDITIONS ENTITY class
  % Constrained entities IDs and values of the BCs are defined by the user.

  properties (Access = public)
    % Boundary condition identifier (mainly for error messages)
    name
    % Total number of constrained source entities
    totSrcEnts = 0
    % Total number of constrained target entities
    totTargetEnts = 0
    % Number of constrained source entities for dof component
    nSourceEnts
    % Raw indices of source constrained entities
    sourceEnts
    % Indices of target constrained entities
    targetEnts
    % Number of constrained target entities for each dof component
    nTargetEnts
    % Number of input times
    nTimes
    % Set of input times
    times
    % struct with data to set the boundary condition in times
    bcData = struct([])
    % Values of currently-stored boundary conditions
    availVals
    % Time id of currently-stored boundary conditions
    availSteps
    % the type of the source entity where the BC is applied
    sourceField
    % the position in space of the source entities
    entityPos
    % logical index to easily deactivate a bc entity
    isActiveEntity
    % source to target geometric mapping operator
    entsMap
  end

  properties (Access = private)
    dof
    readSetFlag = true
  end

  methods (Access = public)
    % Class constructor method
    function obj = BoundaryEntities(name, srcFld)
      % Calling the function to set object properties
      obj.name = name;
      obj.sourceField = srcFld;
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
      obj.bcData = [obj.bcData; struct('time', tVal, 'value', params.value)];

      % reorder time in ascending order
      [obj.times,s] = sort(obj.times,"ascend");
      obj.bcData = obj.bcData(s);

      obj.nTimes = sum(obj.times >= 0.0);

      % check for repeated bc times
      if length(unique(obj.times))~=length(obj.times)
        error("Multiple BC events with same time are" + ...
          "not allowed")
      end

    end


  end

  methods (Access = public)

    function setEntities(obj,type,list,comp,mesh)

      if obj.totSrcEnts > 0
        gresLog().warning(2,['Entities for boundary condition %s ' ...
          'already defined. GReS will overwrite them.'])
      end

      [obj.nSourceEnts, obj.sourceEnts, obj.entityPos] = readEntitySet(obj,type,list,comp,mesh);
      obj.nSourceEnts = reshape(obj.nSourceEnts,1,[]);

      obj.totSrcEnts = sum(obj.nSourceEnts);

      if (obj.totSrcEnts == 0)
        error('No boundary conditions are prescribed for %s BC', obj.name);
      end

      obj.availVals = zeros(obj.totSrcEnts,2);
      obj.availSteps = zeros(2,1);

    end

    function computeTargetEntities(obj,grid,targetField,srcField)

      obj.nTargetEnts = zeros(numel(obj.nSourceEnts),1);

      n = 0;

      comp = find(obj.nSourceEnts > 0);
      comp = reshape(comp,1,[]);

      for i = comp
        % process components individually
        srcID = obj.sourceEnts(n+1:n+obj.nSourceEnts(i));

        n = n+obj.nSourceEnts(i);

        if srcField == entityField.surface && targetField == entityField.cell
          % For FV only: treated differently
          obj.targetEnts = srcID;
          continue
        end

        [inflMap,targEnts] = getIncidenceMap(targetField,grid,srcField,srcID);

        obj.nTargetEnts(i) = numel(targEnts);
        obj.targetEnts = [obj.targetEnts; targEnts];

        % concatenate block diagonal sparse maps for each component set
        obj.entsMap = blkdiag(obj.entsMap,inflMap);

      end

      obj.totTargetEnts = sum(obj.nTargetEnts);

      obj.isActiveEntity = true(obj.totTargetEnts,1);
      % obj.availVals = zeros(obj.totTargetEnts,2);
      % obj.availSteps = zeros(2,1);

    end



    function removeTargetEntities(obj,list)
      % remove BC entities that are contained in an input list
      % ignores entries of list that are not valid entities

      isEntActive = ~ismember(obj.targetEnts,list);
      obj.isActiveEntity(~isEntActive) = false;
      obj.targetEnts(~isEntActive) = [];
      obj.entsMap(~isEntActive,:) = [];

      % update the number of entities
      n = 0;
      ncomp = numel(obj.nTargetEnts);
      l = zeros(ncomp,1);

      for i = 1:ncomp
        l(i) = sum(isEntActive(n+1:obj.nTargetEnts(i)));
        n = n + obj.nTargetEnts(i);
      end

      obj.nTargetEnts = l;
      obj.totTargetEnts = sum(obj.nTargetEnts);
      % obj.availVals = zeros(obj.totTargetEnts,2);

    end

  end



  methods (Access=private)

    function [nEnts, ents, entsPosition] = readEntitySet(obj,type, ents, components, mesh)
      % read entity set and return the number of entities for each
      % component, the list of entities and their reference location

      switch lower(type)
        % input file for list of entities
        case "bclistfile"
          [nEnts, ents] = obj.readListFile(ents);
          entsPosition = getLocation(obj,ents,mesh);
          return
        case {'tags','tag'}
          tags = ents;
          entsID = getEntitiesFromTags(obj.sourceField,mesh,obj.sourceField,tags);
        case "bclist"
          entsID = ents;
        case "box"
          boxSize = ents;
          Lx = boxSize(1:2);
          Ly = boxSize(3:4);
          Lz = boxSize(5:6);
          switch obj.sourceField
            case "node"
              c = mesh.coordinates;
            case "cell"
              c = mesh.cellCentroid;
            otherwise
              error("Error for BC %s: entityListType 'box' is not valid for BC of type %s", obj.name, obj.sourceField)
          end

          entsID = all([ c(:,1) > Lx(1), c(:,1) < Lx(2),...
            c(:,2) > Ly(1), c(:,2) < Ly(2),...
            c(:,3) > Lz(1), c(:,3) < Lz(2)],2);

          entsID = find(entsID);

        otherwise
          error("Unrecognized entityListType '%s' for Boundary condition '%s'.\n" + ...
            "Valid entries for 'entityListType' are: 'bclist','bclistfile','tags','box'.",type,obj.name)
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
      entsID = reshape(entsID,[],1);
      ents = repmat(entsID,sum(compID),1);
      entsPosition = getLocation(obj,ents,mesh);

    end



    function pos = getLocation(obj,ents,mesh)

      switch obj.sourceField
        case "node"
          pos = mesh.coordinates(ents,:);
        case "surface"
          pos = mesh.surfaceCentroid(ents,:);
        case "cell"
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


  methods (Static)

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

  end
end
