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

    vals = getValues(obj,t)

    addBCEvent(obj,varargin)


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

    dataVals = readDataSet(obj,varargin)

    [nEnts, ents, entsPosition] =  ...
      readEntitySet(obj,type, ents, components, mesh)


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

  end


  methods (Static)

    [nEnts,ents] = readListFile(fileName)

  end

end
