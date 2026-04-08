classdef entityField
  % general class for handling field location and interactions

  % possible location of variable fields
  enumeration
    node
    face % not supported yet
    cell
    surface
  end


  methods

    function infl = getIncidenceID(target,grid,source,sourceList)
      % given entities of type 'source', return the connected entities of
      % type 'target'.
      % sourceList - optional input with id of source entities

      % output: list: list of indices connected (with repetitions)
      % ptr: pointer mapping each source id with the first corresponding
      % target entity in the list

      if nargin < 4
        sourceList = [];
      end

      target = entityField(target);

      switch source
        case entityField.node
          infl = getIncidenceIDFromNode(target,grid,sourceList);
        case entityField.surface
          infl = getIncidenceIDFromSurface(target,grid,sourceList);
        case entityField.cell
          infl = getIncidenceIDFromCell(target,grid,sourceList);
      end

    end

    function ents = getEntitiesFromTags(target,grid,source,tags)
      switch source
        case entityField.node
          error("Tag entity retrieval is not valid for source entity of type '%s'",source)
        case entityField.surface
          list = find(ismember(grid.surfaces.tag,tags));
        case entityField.cell
          list = find(ismember(grid.cells.tag,tags));
      end

      ents = getEntitiesList(target,grid,source,list);

    end

    function ents = getEntitiesList(target,grid,varargin)
      if nargin < 3
        source = target;
        list = [];
      elseif nargin < 4
        source = varargin{1};
        list = [];
      else
        source = varargin{1};
        list = varargin{2};
      end

      ents = getIncidenceID(target,grid,source,list);
      ents = unique(ents.getData);
    end

    function nEnts = getNumberOfEntities(target,grid,varargin)
      ents = getEntitiesList(target,grid,varargin{:});
      nEnts = numel(ents);
    end


    function infl = getIncidenceIDFromNode(target,grid,sourceList)

      if isempty(sourceList)
        sourceList = 1:grid.nNodes;
      end

      switch target
        case entityField.node
          infl = ArrayOfArrays(sourceList,ones(numel(sourceList),1));
        otherwise
          error("Incidence from node to %s is not yet available")
      end

    end

    function infl = getIncidenceIDFromSurface(target,grid,sourceList)

      if isempty(sourceList)
        sourceList = 1:grid.surfaces.num;
      end

      switch target
        case entityField.node
          list = grid.getSurfNodes(sourceList);
          rl = grid.surfaces.numVerts(sourceList);
        case entityField.surface
          list = sourceList;
          rl = ones(numel(sourceList),1);
        case entityField.cell
          error("Incidence from surface to %s is not yet available")
      end

      list = list';
      infl = ArrayOfArrays(list(:),rl);

    end

    function infl = getIncidenceIDFromCell(target,grid,sourceList)

      if isempty(sourceList)
        sourceList = 1:grid.cells.num;
      end

      switch target
        case entityField.node
          list = grid.getCellNodes(sourceList);
          rl = grid.cells.numVerts(sourceList);
        case entityField.cell
          list = sourceList;
          rl = ones(numel(sourceList),1);
        otherwise
          error("Incidence from cell to %s is not yet available",target)
      end

      list = list';
      infl = ArrayOfArrays(list(:),rl);

    end


    function [map,varargout] = getIncidenceMap(target,grid,source,srcList)
      % return a sparse matrix  map that perform geometrical interpolation
      % from entity of type 'source' to entity of type 'target'
      % sourceList - optional input with id of source entities

      % output: list: list of indices connected, ptr: pointer mapping each
      % source id with the first corresponding target entity in the list

      if nargin < 4
        srcList = [];
      end


      switch entityField(source)
        case entityField.node
          [map,ents] = getIncidenceMapFromNode(target,grid,srcList);
          % map = getIncidenceMapFromNode(target,grid,srcList);
        case entityField.surface
          [map,ents] = getIncidenceMapFromSurface(target,grid,srcList);
        case entityField.cell
          [map,ents] = getIncidenceMapFromCell(target,grid,srcList);
      end

      if nargout > 1
        varargout{1} = ents;
      end

    end

    function [map,targEnts] = getIncidenceMapFromNode(target,grid,srcList)

      if isempty(srcList)
        srcList = 1:grid.nNodes;
      end

      switch target
        case entityField.node
          n = numel(srcList);
          map = speye(n);
          targEnts = srcList;
        otherwise
          error("Incidence map from node to %s is not yet available")
      end

    end


    function [map,targEnts] = getIncidenceMapFromSurface(target,grid,srcList)

      if isempty(srcList)
        srcList = 1:grid.surfaces.num;
      end

      % need to update logic for node/volume influence

      switch target
        case entityField.node
          [nodes,areas] = grid.getNodeInfluence(entityField.surface,srcList); % nodes is an ArrayOfArrays
          [nList,ptr] = getData(nodes);   
          [targEnts,~,ii] = unique(nList); 
          jj = repelem((1:numel(srcList))',diff(ptr),1);
          map = sparse(ii,jj,areas,numel(targEnts),numel(srcList));

        case entityField.surface
          n = numel(srcList);
          map = speye(n);
          targEnts = srcList;

        case entityField.cell
          % this combination is only needed by FV solvers where influence
          % map is not used
          map = [];
          targEnts = [];
          %error("Incidence map from surface to %s is not yet available")
      end

    end


    function [map,targEnts] = getIncidenceMapFromCell(target,grid,srcList)

      if isempty(srcList)
        srcList = 1:grid.cells.num;
      end

      switch target
        case entityField.node

          [nodes,vols] = grid.getNodeInfluence(entityField.cell,srcList); % nodes is an ArrayOfArrays
          [nList,ptr] = getData(nodes);
          [targEnts,~,ii] = unique(nList);
          jj = repelem((1:numel(srcList))',diff(ptr),1);
          map = sparse(ii,jj,vols,numel(targEnts),numel(srcList));

        case entityField.cell
          n = numel(srcList);
          map = speye(n);
          targEnts = srcList;
        otherwise
          error("Incidence map from cell to %s is not yet available")
      end
    end


    function size = getEntitySize(ent,msh,list)

      if nargin < 3
        list = getIncidenceID(ent,msh,ent);
      end

      switch ent
        case entityField.node
          size = ones(numel(list),1);
        case entityField.surface
          size = msh.surfaces.area(list);
        case entityField.cell
          size = msh.cells.volume(list);
      end
    end

  end

end



