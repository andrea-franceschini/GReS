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

      if strcmp(obj.grid.topology.meshType,"Unstructured")
        % Calling the function to read input data from file
        obj.addBCs(varargin{2:end});
        obj.computeBoundaryProperties(obj.grid);
      end
    end

    function addBCs(obj,varargin)

      % add a list of multiple boundary conditions from file
      assert(isscalar(varargin),"addBCs method is valid only with scalar" + ...
        " xml file name or ")
      input = readInput(varargin{:});

      for i = 1:numel(input.BC)
        addBC(obj,input.BC(i));
      end

    end

    function addBC(obj,varargin)

      default = struct('name',[],...
        'variable',[],...
        'targetEntity',[],...
        'entityList',[],...
        'entityListType',[],...
        'type',[],...
        'components',missing);

      params = readInput(default,varargin{:});

      name = params.name;
      if obj.db.isKey(name)
        error("'%s' boundary condition name already defined", name);
      end

      targetEnt = lower(params.targetEntity);
      variable = params.variable;
      type = params.type;


      if (~ismember(string(targetEnt), ["node", "surface", "cell", "volumeforce"]))
        error(['%s condition is unknown\n', ...
          'Accepted types are: node   -> Boundary cond. on nodes\n',...
          '                    surface   -> Boundary cond. on surfaces\n',...
          '                    cell   -> Boundary cond. on elements\n',...
          '                    volumeforce -> Volume force on elements'], targetEnt);
      end


      if ~strcmp(targetEnt,"VolumeForce")
        if (~ismember(type, ["Dirichlet", "Neumann", "Seepage"]))
          error(['Error in BC %s : %s boundary condition is not admitted\n', ...
            'Accepted types are: Dirichlet, Neumann, Seepage'], name, type);
        end
      end


      switch targetEnt
        case 'volumeforce'
          bc = struct('data', [], ...
            'cond',targetEnt, 'variable', variable);
        case {'surface','node','cell'}
          bc = struct('data', [], ...
            'cond',targetEnt,'type', type, 'variable', variable);
        otherwise
          error('Unrecognized BC item %s for Boundary condition %s', ...
            targetEnt, name)
      end

      if ismissing(params.components)
        comp = [];
      else
        comp = params.components;
      end

      bcEnt = BoundaryEntities(name,targetEnt);
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

    function vals = getVals(obj, identifier, t)
      vals = obj.getData(identifier).data.getValues(t);
    end

    function cond = getCond(obj, identifier)
      cond = obj.getData(identifier).cond;
    end

    function name = getName(obj, identifier)
      name = obj.getData(identifier).data.name;
    end

    function type = getType(obj, identifier)
      if ~strcmp(obj.getCond(identifier),'VolumeForce')
        type = obj.getData(identifier).type;
      else
        type = 'VolumeForce';
      end
    end

    function variable = getVariable(obj, identifier)
      variable = obj.getData(identifier).variable;
    end

    function dofs = getCompEntities(obj,identifier,ents)
      % get component dof of Dirichlet BC loaded entities
      nEnts = getNumbLoadedEntities(obj,identifier);
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

    function ents = getEntities(obj,identifier)
      ents = obj.getData(identifier).data.entities;
    end

    function ents = getBCentities(obj,identifier)
      if isfield(obj.getData(identifier),"loadedEnts")
        ents = getLoadedEntities(obj,identifier);
      else
        ents = getEntities(obj,identifier);
      end
    end

    function ents = getLoadedEntities(obj, identifier)
      % return loaded entities for Volume or Surface BCs
      ents = obj.getData(identifier).loadedEnts;
    end

    function nEnts = getNumbEntities(obj,identifier)
      nEnts = obj.getData(identifier).data.nEntities;
    end

    function ents = getNumbLoadedEntities(obj, identifier)
      bc = obj.getData(identifier);
      if isfield(bc,'nloadedEnts')
        % Surface BC
        ents = bc.nloadedEnts;
      else
        % Node BC
        ents = bc.data.nEntities;
      end
    end

    function infl = getEntitiesInfluence(obj, identifier)
      infl = [];
      if isfield(obj.getData(identifier),"entitiesInfl")
        infl = obj.getData(identifier).entitiesInfl;
      end
    end

    function setDofs(obj, identifier, list)
      obj.getData(identifier).data.entities = list;
    end

    function computeBoundaryProperties(obj, grid)

      % preprocess surface/volume boundary conditions

      msh = grid.topology;
      elem = grid.cells;

      keys = obj.db.keys;

      for bcId = 1:length(keys)
        key = keys{bcId};
        cond = obj.getCond(key);
        phys = obj.getVariable(key);

        if any(strcmp(cond, ["VolumeForce","SurfBC"]))

          ents = obj.getEntities(key);
          nEnts = obj.getData(key).data.nEntities;
          nLoadEnts = zeros(numel(nEnts),1);
          loadedEnts = [];
          entsInfl = [];

          N = 0;
          for i = 1:numel(nEnts)
            ents_i = ents(N+1:N+nEnts(i));
            if strcmp(cond,'VolumeForce')
              tmpMat = msh.cells(ents_i, :)';
              nEntries = sum(msh.cellNumVerts(ents_i));
            else
              tmpMat = msh.surfaces(ents_i, :)';
              nEntries = sum(msh.surfaceNumVerts(ents_i));
            end

            loadedEnts_i = unique(tmpMat(tmpMat ~= 0));
            nLoadEnts(i) = numel(loadedEnts_i);

            % Preallocate row,col,val indices for sparse assembly
            %           n = sum(msh.cellNumVerts())
            [r,c,v] = deal(zeros(nEntries,1));
            k = 0;
            for j = 1:nEnts(i)
              el = ents_i(j);
              if strcmp(cond,'VolumeForce')
                nodInf = findNodeVolume(elem,el);
                nodes = msh.cells(el,:);
              else
                nodInf = findNodeArea(elem,el);
                nodes = msh.surfaces(el,:);
              end
              loadEntsLoc = find(ismember(loadedEnts_i,nodes));
              nn = numel(nodInf);
              r(k+1:k+nn) = loadEntsLoc;
              c(k+1:k+nn) = repelem(j,nn);
              v(k+1:k+nn) = nodInf;
              k = k + nn;
            end
            entsInfl = blkdiag(entsInfl,sparse(r,c,v));
            N = N + nEnts(i);
            loadedEnts = [loadedEnts; loadedEnts_i];
          end

          if strcmp(obj.getType(key), 'Dirichlet')
            entsInfl = entsInfl./sum(entsInfl,2);
          end

          % update bc struct with additional properties
          entry = obj.getData(key);
          entry.entitiesInfl = entsInfl;
          entry.loadedEnts = loadedEnts;
          entry.nloadedEnts = nLoadEnts;
          obj.db(key) = entry;
        end
      end

    end


    function removeBCentities(obj,bcId,list)
      % remove BC entities that are contained in an input list
      % ignores entries of list that are not valid entities

      bc = getData(obj,bcId);

      if isfield(bc,"loadedEnts")
        % remove loaded entity
        loadedEnts = getLoadedEntities(obj,bcId);
        isEntActive = ~ismember(loadedEnts,list);
        bc.loadedEnts = bc.loadedEnts(isEntActive);

        % remove corresponding rows of entity influence matrix
        bc.entitiesInfl(~isEntActive,:) = [];

        % update the number of loaded entities
        n = 0;
        ncomp = numel(bc.nloadedEnts);
        l = zeros(ncomp,1);
        for i = 1:ncomp
          l(i) = sum(isEntActive(n+1:bc.nloadedEnts(i)));
          n = n + bc.nloadedEnts(i);
        end

        bc.nloadedEnts = l;

      else
        % remove entities
        bcEnts = bc.data;
        isEntActive = ~ismember(bcEnts.entities,list);
        bcEnts.isActiveEntity(~isEntActive) = false;
        bcEnts.entities(~isEntActive) = [];

        % update the number of entities
        n = 0;
        ncomp = numel(bcEnts.nEntities);
        l = zeros(ncomp,1);

        for i = 1:ncomp
          l(i) = sum(isEntActive(n+1:bcEnts.nEntities(i)));
          n = n + bcEnts.nEntities(i);
        end

        bcEnts.nEntities = l;
        bcEnts.totEnts = sum(bcEnts.nEntities);
        bcEnts.availVals = zeros(bcEnts.totEnts,2);

      end

      obj.db(bcId) = bc;

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