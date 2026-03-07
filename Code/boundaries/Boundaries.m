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

    function addBCs(obj,varargin)

      % add a list of multiple boundary conditions from file
      assert(isscalar(varargin),"addBCs method is valid only with scalar" + ...
        " xml file name or struct")
      input = readInput(varargin{:});

      for i = 1:numel(input.BC)
        addBC(obj,input.BC(i));
      end

    end

    function addBC(obj,varargin)

      default = struct('name',string.empty,...
        'variable',string.empty,...
        'field',string.empty,...
        'entityList',double.empty,...
        'entityListType',string.empty,...
        'type',string.empty,...
        'components',missing,...
        'essential',missing);

      params = readInput(default,varargin{:});

      name = params.name;
      if obj.db.isKey(name)
        error("'%s' boundary condition name already defined", name);
      end

      srcEnt = entityField(lower(params.field));
      variable = lower(params.variable);
      type = lower(params.type);

      % BC type
      % THE USER IS RESPONSIBLE FOR INPUTTING THE CORRECT BCTYPE NAME.
      % THIS IS ONLY CHECKED AT THE PHYSICS_SOLVER LEVEL

      essentialFlag = false;

      if BCtype.isCustomBC(type)
        if ~ismissing(params.essential)
          essentialFlag = logical(params.essential);
        else
          gresLog().warning(1, sprintf( ...
            "Custom boundary condition '%s' has been automatically considered to be NOT essential!\n" + ...
            "To override this behavior, specify the logical field 'essential' in the Boundary condition input.", ...
            bctype));
        end
      elseif BCtype(type) == BCtype.dirichlet
        essentialFlag = true;
      end

      bc = struct('data', [], 'essential',essentialFlag,...
        'sourceField',srcEnt,'type', type, 'variable', variable);

      if ismissing(params.components)
        comp = [];
      else
        comp = params.components;
      end

      bcEnt = BoundaryEntities(name,srcEnt);
      bcEnt.setEntities(params.entityListType,...
        params.entityList,...
        comp,obj.grid.topology);

      bc.data = bcEnt;

      % add BC to the database
      obj.db(name) = bc;

      % add the bc event
      if isfield(params,"BCevent")
        for i = 1:numel(params.BCevent)
          addBCEvent(obj,name,params.BCevent(i))
        end
      end

      setBCList(obj);

      % finalize the boundary condition
      % obj.computeBoundaryProperties(name);

    end


    function addBCEvent(obj,bcId,varargin)

      obj.getData(bcId).data.addBCEvent(varargin{:});

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
      ents = obj.getData(identifier).data.sourceEntities;
    end

    function ents = getTargetEntities(obj,identifier)
      % get raw entities as specified in the BC input
      ents = obj.getData(identifier).data.targetEntities;
    end

    function dofs = getBCdofs(obj,bcId)
      % get the constrained degree-of-freedom of the specified bcs

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
      bcEnt = obj.getData(bcId).data;
      bcEnt.computeTargetEntities(obj.grid,target,src);

    end


    function nEnts = getNumbSourceEntities(obj,identifier)
      nEnts = obj.getData(identifier).data.nSrcEntities;
    end

    function nEnts = getNumbTargetEntities(obj,identifier)
      nEnts = obj.getData(identifier).data.nTargetEntities;
    end

    function infl = getEntitiesInfluence(obj, identifier)
      infl = obj.getData(identifier).data.entsMap;
    end

    % function setDofs(obj, identifier, list)
    %   getData(obj,identifier).data.sourceEntities = list;
    % end

    %   function computeBoundaryProperties(obj,bcId)
    %
    %     % preprocess boundary conditions once the target entity is known
    %
    %     msh = obj.grid.topology;
    %     elem = obj.grid.cells;
    %
    %     % target entity field
    %     source = obj.getBCfield(bcId);
    %
    %       ents = obj.getEntities(bcId);
    %       nSrcEnts = obj.getData(bcId).data.nEntities;
    %       nTargetEnts = zeros(numel(nSrcEnts),1);
    %
    %       loadedEnts = [];
    %       entsInfl = [];
    %
    %       N = 0;
    %
    %       for i = 1:numel(nSrcEnts)
    %         % process components individually
    %         ents_i = ents(N+1:N+nSrcEnts(i));
    %         inflMap_i = getEntitiesInterpolation()
    %         if strcmp(source,'volumeforce')
    %           tmpMat = msh.cells(ents_i, :)';
    %           nEntries = sum(msh.cellNumVerts(ents_i));
    %         else
    %           tmpMat = msh.surfaces(ents_i, :)';
    %           nEntries = sum(msh.surfaceNumVerts(ents_i));
    %         end
    %
    %         loadedEnts_i = unique(tmpMat(tmpMat ~= 0));
    %         nTargetEnts(i) = numel(loadedEnts_i);
    %
    %         % Preallocate row,col,val indices for sparse assembly
    %         %           n = sum(msh.cellNumVerts())
    %         [r,c,v] = deal(zeros(nEntries,1));
    %         k = 0;
    %         for j = 1:nSrcEnts(i)
    %           el = ents_i(j);
    %           if strcmp(source,'volumeforce')
    %             nodInf = findNodeVolume(elem,el);
    %             nodes = msh.cells(el,:);
    %           else
    %             nodInf = findNodeArea(elem,el);
    %             nodes = msh.surfaces(el,:);
    %           end
    %           loadEntsLoc = find(ismember(loadedEnts_i,nodes));
    %           nn = numel(nodInf);
    %           r(k+1:k+nn) = loadEntsLoc;
    %           c(k+1:k+nn) = repelem(j,nn);
    %           v(k+1:k+nn) = nodInf;
    %           k = k + nn;
    %         end
    %         entsInfl = blkdiag(entsInfl,sparse(r,c,v));
    %         N = N + nSrcEnts(i);
    %         loadedEnts = [loadedEnts; loadedEnts_i];
    %       end
    %
    %       if strcmpi(obj.getType(bcId), 'dirichlet')
    %         entsInfl = entsInfl./sum(entsInfl,2);
    %       end
    %
    %       % update bc struct with additional properties
    %       entry = obj.getData(bcId);
    %       entry.entitiesInfl = entsInfl;
    %       entry.targetEntities = loadedEnts;
    %       entry.nTargetEnts = nTargetEnts;
    %       obj.db(bcId) = entry;
    %     end
    % end



    function removeBCentities(obj,bcId,list)
      % remove BC entities that are contained in an input list
      % ignores entries of list that are not valid entities

      getData(obj,bcId).removeTargetEntities(list);

    end

    function out = isEssential(obj,bcId)

      out = obj.getData(bcId).essential;

    end



    function bcList = getBCList(obj)
      % Get list of boundary condition names ensuring that Dirichlet are
      % applied at last

      if isempty(obj.bcList)
        setBCList(obj);
      end

      bcList = obj.bcList;
    end
  end

  methods (Access = private)

    function setBCList(obj)
      % set correct order of boundary conditions. Dirichlet last
      bcTypeList = [];
      for bcId = string(obj.db.keys)
        bcTypeList = [bcTypeList, obj.getType(bcId)];
      end
      idxDir  = strcmp(bcTypeList, 'Dirichlet');
      bcOrd = [find(~idxDir), find(idxDir)];
      bcNames = obj.db.keys;
      obj.bcList = string(bcNames(bcOrd));
    end

  end

end