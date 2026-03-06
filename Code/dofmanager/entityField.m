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

      switch source
        case entityField.node
          [list,ptr] = getIncidenceIDFromNode(target,mesh,sourceList);
        case entitiyField.surface
          [list,ptr] = getIndicenceIDFromSurface(target,mesh,sourceList);
        case entityField.cell
          [list,ptr] = getIndicenceIDFromCell(target,mesh,sourceList);
      end

      varargout{1} = list;

      if nargout > 1
        varargout{2} = ptr;
      end
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

    function [list,ptr] = getIndicenceIDFromSurface(target,mesh,sourceList)

      if isempty(sourceList)
        sourceList = 1:mesh.nSurfaces;
      end

      switch target
        case entityField.node
          ptr = mesh.surfaceNumVerts(sourceList);
          tmp = mesh.surfaces(sourceList,:);
          list = reshape(tmp',[],1);
        case entityField.surface
          list = sourceList;
          ptr = reshape(1:numel(list),[],1);
        case entityField.cell
          error("Incidence from surface to %s is not yet available")
      end

    end

    function [list,ptr] = getIndicenceIDFromCell(target,mesh,sourceList)

      if isempty(sourceList)
        sourceList = 1:mesh.nCells;
      end

      switch target
        case entityField.node
          ptr = [1; cumsum(mesh.cellNumVerts(sourceList))];
          tmp = mesh.cells(sourceList,:);
          list = reshape(tmp',[],1);
        case entityField.cell
          list = sourceList;
          ptr = reshape(1:numel(list),[],1);
        otherwise
          error("Incidence from cell to %s is not yet available",target)
      end

    end


    function map = getIncidenceMap(target,grid,source,srcList)
      % return a sparse matrix  map that perform geometrical interpolation
      % from entity of type 'source' to entity of type 'target'
      % sourceList - optional input with id of source entities

      % output: list: list of indices connected, ptr: pointer mapping each
      % source id with the first corresponding target entity in the list

      if nargin < 4
        srcList = [];
      end


      switch source
        case entityField.node
          error("Incidence map from node is not yet available")
          % map = getIncidenceMapFromNode(target,grid,srcList);
        case entitiyField.surface
          map = getIncidenceMapFromSurface(target,grid,srcList);
        case entityField.cell
          map = getIncidenceMapFromCell(target,grid,srcList);
      end

    end


    function [map,varargout] = getIncidenceMapFromSurface(target,grid,srcList)
      switch target
        case entityField.node
          [nodeID,ptr] = getIncidenceIDFromSurface(target,grid,srcList);
          [targEnts,~,n] = unique(nodeID);
          m = assembler(numel(nodeID),numel(targEnts),numel(srcList));
          k = 0;
          for id = srcList'
            k = k+1;
            A = findNodeArea(grid.cell,id);
            ii = n(ptr(k:k+1));
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

      if nargout > 1
        varargout{1} = targEnts;
      end

    end


    function map = getIncidenceMapFromCell(target,grid,srcList)
      switch target
        case entityField.node
          [nodeID,ptr] = getIncidenceIDFromCell(target,grid,srcList);
          [nodes,~,n] = unique(nodeID);
          m = assembler(numel(nodeID),numel(nodes),numel(srcList));
          k = 0;
          for id = srcList'
            k = k+1;
            V = findNodeVolume(grid.cell,id);
            ii = n(ptr(k:k+1));
            m.localAssembly(ii,k,V./sum(V));
          end
          map = m.sparseAssembly();
        case entityField.cell
          n = numel(srcList);
          map = speye(n);
        otherwise
          error("Incidence map from cell to %s is not yet available")
      end
    end

  end

end



