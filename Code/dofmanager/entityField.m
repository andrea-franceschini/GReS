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

    function varargout = getIncidenceID(target,mesh,source,sourceList)
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
          [list,ptr] = getIncidenceIDFromNode(target,mesh,sourceList);
        case entityField.surface
          [list,ptr] = getIncidenceIDFromSurface(target,mesh,sourceList);
        case entityField.cell
          [list,ptr] = getIncidenceIDFromCell(target,mesh,sourceList);
      end

      varargout{1} = list;

      if nargout > 1
        varargout{2} = ptr;
      end
    end

    function ents = getEntitiesFromTags(target,mesh,source,tags)
      switch source
        case entityField.node
          error("Tag entity retrieval is not valid for source entity of type '%s'",source)
        case entityField.surface
          list = find(ismember(mesh.surfaceTag,tags));
        case entityField.cell
          list = find(ismember(mesh.cellTag,tags));
      end

      ents = getEntitiesList(target,mesh,source,list);

    end

    function ents = getEntitiesList(target,mesh,varargin)
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

      ents = getIncidenceID(target,mesh,source,list);
      ents = unique(ents);
    end

    function nEnts = getNumberOfEntities(target,mesh,varargin)
      ents = getEntitiesList(target,mesh,varargin{:});
      nEnts = numel(ents);
    end


    function [list,ptr] = getIncidenceIDFromNode(target,mesh,sourceList)

      if isempty(sourceList)
        sourceList = 1:mesh.nNodes;
      end

      switch target
        case entityField.node
          list = sourceList;
          ptr = reshape(1:numel(list),[],1);
        otherwise
          error("Incidence from node to %s is not yet available")
      end

    end

    function [list,ptr] = getIncidenceIDFromSurface(target,mesh,sourceList)

      if isempty(sourceList)
        sourceList = 1:mesh.nSurfaces;
      end

      switch target
        case entityField.node
          ptr = [0;cumsum(mesh.surfaceNumVerts(sourceList))];
          tmp = mesh.surfaces(sourceList,:);
          list = reshape(tmp',[],1);
        case entityField.surface
          list = sourceList;
          ptr = reshape(1:numel(list),[],1);
        case entityField.cell
          error("Incidence from surface to %s is not yet available")
      end

    end

    function [list,ptr] = getIncidenceIDFromCell(target,mesh,sourceList)

      if isempty(sourceList)
        sourceList = 1:mesh.nCells;
      end

      switch target
        case entityField.node
          ptr = [0; cumsum(mesh.cellNumVerts(sourceList))];
          tmp = mesh.cells(sourceList,:);
          list = reshape(tmp',[],1);
        case entityField.cell
          list = sourceList;
          ptr = reshape(1:numel(list),[],1);
        otherwise
          error("Incidence from cell to %s is not yet available",target)
      end

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

      target = entityField(target);


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
        srcList = 1:grid.topology.nNodes;
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
        srcList = 1:grid.topology.nSurfaces;
      end

      switch target
        case entityField.node
          [nodeID,ptr] = getIncidenceIDFromSurface(target,grid.topology,srcList);
          [targEnts,~,n] = unique(nodeID);
          m = assembler(numel(nodeID),numel(targEnts),numel(srcList));
          k = 0;
          for id = srcList'
            k = k+1;
            A = findNodeArea(grid.cells,id);
            ii = n(ptr(k)+1:ptr(k+1));
            m.localAssembly(ii,k,A./sum(A));
          end
          map = m.sparseAssembly();
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
        srcList = 1:grid.topology.nCells;
      end

      switch target
        case entityField.node
          [nodeID,ptr] = getIncidenceIDFromCell(target,grid.topology,srcList);
          [targEnts,~,n] = unique(nodeID);
          m = assembler(numel(nodeID),numel(targEnts),numel(srcList));
          k = 0;
          for id = srcList'
            k = k+1;
            V = findNodeVolume(grid.cells,id);
            ii = n(ptr(k)+1:ptr(k+1));
            m.localAssembly(ii,k,V./sum(V));
          end
          map = m.sparseAssembly();
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
          size = ones(list,1);
        case entityField.surface
          size = msh.surfaceArea(list);
        case entityField.cell
          size = msh.cellVolume(list);
      end
    end

  end

end



