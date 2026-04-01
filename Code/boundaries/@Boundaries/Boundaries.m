classdef Boundaries < handle
  % BOUNDARY CONDITIONS - General boundary conditions class

  properties (Access = public)
    % Creation of a Map object for the boundary conditions
    db
  end

  properties (Access = private)
    grid
    bcList
    % to ensure neumann is applied before dirichlet in case of overlapping
    % entity definition
  end

  methods (Access = public)


    function obj = Boundaries(varargin)
      % MATLAB evaluates the assignment expression for each instance, which
      % ensures that each instance has a unique value
      obj.db = containers.Map('KeyType','char','ValueType','any');

      if nargin == 0
        return
      end

      if nargin == 1
        obj.grid = varargin{1};
        return
      end

      obj.grid = varargin{1};

      % % check grid
      % assert()

      if strcmp(obj.grid.topology.meshType,"Unstructured")
        % Calling the function to read input data from file
        obj.addBCs(varargin{2:end});
      end
    end


    addBC(obj,varargin)

    addBCEvent(obj,varargin)


    function addBCs(obj,varargin)

      % add a list of multiple boundary conditions from file
      assert(isscalar(varargin),"addBCs method is valid only with scalar" + ...
        " xml file name or struct")
      input = readInput(varargin{:});

      for i = 1:numel(input.BC)
        addBC(obj,input.BC(i));
      end

    end



    function bc = getData(obj,identifier)
      % Check if the identifier defined by the user is a key of the Map object
      if (obj.db.isKey(identifier))
        bc = obj.db(identifier);
      else
        % Displaying error message if the identifier does not refer
        % to any existing class
        error('Boundary condition %s not present', identifier);
      end
    end

    function vals = getVals(obj, bcId, t)
      % get value of the boundary condition at the target field
      
      % geometric incidence map from source to target bc entities
      M = obj.getEntitiesInfluence(bcId);

      valSrc = getSourceVals(obj,bcId,t);

      type = obj.getType(bcId);

      msg = "The base getVals() method cannot be used for custom boundary condition '" + type + "'." + newline + ...
        "It must be overridden in a method of the PhysicsSolver using it" + newline + ...
        "getVals() only applies to BC types: " + ...
        strjoin(string(enumeration('BCtype')), ", ");

      assert(~BCtype.isCustomBC(type), msg);

      if isEssential(obj,bcId)
          % make map interpolative
          M = M./sum(M,2);
      end

      if type == BCtype.source || type == BCtype.neumann
        % multiply neumann and source bcs by the size of the source entity
        ents = obj.getSourceEntities(bcId);
        S = getEntitySize(obj.getField(bcId),obj.grid.topology,ents);
        valSrc = S.*valSrc;
      end

      % map values from source entity to target entity
      vals = M * valSrc;
    end

    function vals = getSourceVals(obj,bcId,t)
        vals = obj.getData(bcId).data.getValues(t);
    end


    function cond = getField(obj, identifier)
      cond = obj.getData(identifier).sourceField;
    end

    function name = getName(obj, identifier)
      name = obj.getData(identifier).data.name;
    end

    function type = getType(obj, identifier)
      type = obj.getData(identifier).type;
    end

    function variable = getVariable(obj, identifier)
      variable = obj.getData(identifier).variable;
    end

    function dofs = getCompEntities(obj,identifier,ents)
      % get component dof of Dirichlet BC loaded entities
      nEnts = getNumbTargetEntities(obj,identifier);
      % component multiplication of BC entities
      dim = length(nEnts);
      i1 = 1;
      dofs = zeros(numel(ents),1);
      for i = 1 : dim
        i2 = i1 + nEnts(i);
        dofs(i1:i2-1) = dim*(ents(i1:i2-1)-1) + i;
        i1 = i2;
      end
    end

    function ents = getSourceEntities(obj,identifier)
      % get raw entities as specified in the BC input
      ents = obj.getData(identifier).data.sourceEnts;
    end

    function ents = getTargetEntities(obj,identifier)
      % get raw entities as specified in the BC input
      ents = obj.getData(identifier).data.targetEnts;
    end

    function dofs = getDofs(obj,bcId,dofm)
      % get the constrained degree-of-freedom of the specified bcs

      nEnts = getNumbTargetEntities(obj,bcId);
      ents = getTargetEntities(obj,bcId);

      % transform entities in dof numbering
      if nargin > 2
        var = obj.getVariable(bcId);
        varId = dofm.getVariableId(var);
        ents = getLocalEnts(dofm,varId,ents);
      end

      % component multiplication
      dim = length(nEnts);
      i1 = 1;
      dofs = zeros(numel(ents),1);
      for i = 1 : dim
        i2 = i1 + nEnts(i);
        dofs(i1:i2-1) = dim*(ents(i1:i2-1)-1) + i;
        i1 = i2;
      end
    end

    function dofs = getStateDofs(obj,bcId)
      % get the id of degree-of-freedom of the specified bcs in the state
      % array

      nEnts = getNumbTargetEntities(obj,bcId);
      ents = getTargetEntities(obj,bcId);

      % component multiplication of BC entities
      dim = length(nEnts);
      i1 = 1;
      dofs = zeros(numel(ents),1);
      for i = 1 : dim
        i2 = i1 + nEnts(i);
        dofs(i1:i2-1) = dim*(ents(i1:i2-1)-1) + i;
        i1 = i2;
      end
    end

    function computeTargetEntities(obj,bcId,target)
      % finalize the boundary conditions specifying the entity where the bc
      % is applied

      src = obj.getField(bcId);
      bc = obj.getData(bcId);
      bc.data.computeTargetEntities(obj.grid,target,src);

    end


    function nEnts = getNumbSourceEntities(obj,identifier)
      nEnts = obj.getData(identifier).data.nSorceEnts;
    end

    function nEnts = getNumbTargetEntities(obj,identifier)
      nEnts = obj.getData(identifier).data.nTargetEnts;
    end

    function infl = getEntitiesInfluence(obj, identifier)
      infl = obj.getData(identifier).data.entsMap;
    end

    function setEntities(obj, identifier, list)
      obj.getData(identifier).data.sourceEnts = list;
    end


    function removeTargetEntities(obj,bcId,list)
      % remove BC entities that are contained in an input list
      % ignores entries of list that are not valid entities

      getData(obj,bcId).data.removeTargetEntities(list);

    end

    function out = isEssential(obj,bcId)

      out = obj.getData(bcId).essential;

    end

    function setBCList(obj)
      % set correct order of boundary conditions. Dirichlet last
      bcNames = string(obj.db.keys);
      bcTypeList = strings(numel(bcNames),1);
      for i = 1:numel(bcNames)
        bcId = bcNames(i);
        bcTypeList(i) = string(obj.getType(bcId));
      end
      idxDir  = strcmpi(bcTypeList, 'Dirichlet');
      bcOrd = [find(~idxDir); find(idxDir)];
      bcNames = obj.db.keys;
      obj.bcList = string(bcNames(bcOrd));
    end



    function bcList = getBCList(obj)
      % Get list of boundary condition names ensuring that Dirichlet are
      % applied at last

      bcList = obj.bcList;
    end

    function scaleBC(obj,identifier,scalingFactor)
      % Apply a scaling factor

      obj.getData(identifier).data.bcScale = scalingFactor;
    end
  end

  methods (Access = private)


  end

end
