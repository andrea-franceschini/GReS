classdef InterfaceGrid < Grid
  % Grid class for interfaces

  properties (SetAccess = private)
    neighbors = zeros(0,3)       % [interface cell, slave cell, master cell]
    traceGrids = [Grid(),Grid()]
    intersections = []
  end

  properties (GetAccess = public, SetAccess = private)
    isConforming = false(1,2)
  end

  methods

    function obj = InterfaceGrid(gInterf,g1,g2,varargin)

      obj.setGrid(gInterf);

      if nargin == 1 && isa(gInterf,'InterfaceGrid')
        % Copy the candidate triplets without repeating the contact search.
        obj.traceGrids = gInterf.traceGrids;
        obj.neighbors = gInterf.neighbors;
        obj.intersections = gInterf.intersections;
        obj.isConforming = gInterf.isConforming;
        return
      end
      obj.traceGrids = [Grid();Grid()];
      obj.traceGrids(1).setGrid(g1);
      obj.traceGrids(2).setGrid(g2);

      assert(obj.nDim,"Interface grid must be lower dimensional")

      [obj.neighbors,obj.isConforming] = obj.getNeighbors(varargin{:});

    end


    function [neighbors,conforming] = getNeighbors(obj,varargin)
      % Join the two contact searches on their shared interface cell.

      slaveSearch = ContactSearching(obj,obj.traceGrids(1),varargin{:});
      masterSearch = ContactSearching(obj,obj.traceGrids(2),varargin{:});

      conforming = [slaveSearch.isConforming,masterSearch.isConforming];

      slave = slaveSearch.getConnectivity();
      master = masterSearch.getConnectivity();

      % store triplet of connected facets
      if all(conforming)
        [slaveIds,cells] = find(slave.');
        [masterIds,~] = find(master.');
        neighbors = [cells(:),slaveIds(:),masterIds(:)];
      elseif conforming(1)
        [cells,masterIds] = find(master);
        [slaveIds,~] = find(slave.');
        neighbors = sortrows([cells(:),reshape(slaveIds(cells),[],1),masterIds(:)]);
      elseif conforming(2)
        [cells,slaveIds] = find(slave);
        [masterIds,~] = find(master.');
        neighbors = sortrows([cells(:),slaveIds(:),reshape(masterIds(cells),[],1)]);
      else
        slave = slave.';
        master = master.';
        active = find(any(slave,1) & any(master,1));
        triplets = cell(numel(active),1);
        for k = 1:numel(active)
          cellId = active(k);
          [m,s] = ndgrid(find(master(:,cellId)),find(slave(:,cellId)));
          triplets{k} = [repmat(cellId,numel(s),1),s(:),m(:)];
        end
        neighbors = zeros(0,3);
        if ~isempty(triplets)
          neighbors = vertcat(triplets{:});
        end
      end


    end


    function removeNeighbors(obj,id)

      % Remove triplets selected by a logical mask or a list of row IDs.
      n = size(obj.neighbors,1);

      if islogical(id)
        assert(numel(id) == n && (isvector(id) || isempty(id)), ...
          'InterfaceGrid:InvalidNeighborMask','Expected one logical value per triplet.');
        keep = ~id(:);
      else
        keep = true(n,1);
        if ~isempty(id)
          validateattributes(id,{'numeric'},{'vector','integer','>=',1,'<=',n});
          keep(id) = false;
        end
      end
      newNeighs = obj.neighbors(keep,:);

      if isempty(newNeighs)
        error('No remaining connections in interface grid')
      end

      grids = {obj,obj.traceGrids(1),obj.traceGrids(2)};
      
      for side = 1:3
        grid = grids{side};
        used = false(grid.cells.num,1);
        used(newNeighs(:,side)) = true;

        newIds = cumsum(used);
        newNeighs(:,side) = newIds(newNeighs(:,side));

        if ~all(used)
          grids{side} = grid.getCellGrid(used);
        end
      end

      % update intersection if already computed
      intersect = obj.intersections;
      if ~isempty(intersect) && ~all(keep)
        intersect.center = intersect.center(keep,:);
        intersect.area = intersect.area(keep);
        intersect.verts = intersect.verts.getRows(find(keep));
      end

      obj.setGrid(grids{1});
      obj.traceGrids = [grids{2};grids{3}];
      obj.neighbors = newNeighs;
      obj.intersections = intersect;
    end


    function intersections = findGridIntersections(obj)

      % Cache polygons in candidate triplet order. Quadrature prunes the
      % triplets once their weights are known.
      if ~isempty(obj.intersections)
        intersections = obj.intersections;
        return
      end

      n = size(obj.neighbors,1);
      dim = size(obj.coordinates,2);
      clipSides = find(~obj.isConforming);
      centers = nan(n,dim);
      areas = zeros(n,1);
      counts = zeros(n,1);
      polygons = cell(n,1);

      for k = 1:n 
        % process each triplet

        ids = obj.neighbors(k,:);

        % get vertices of interface face in 2D
        verts = InterfaceGrid.cellPolygon(obj,ids(1));
        project = ~isempty(clipSides) && size(verts,2) == 3;
        if project
          origin = verts(1,:);
          normal = computePolygonNormal(verts);
          [verts,d1,d2] = pointToSurfaceProjection(origin,normal,verts);
          basis = [d1,d2];
        end

        for side = clipSides

          % project side facet onto interface
          clip = InterfaceGrid.cellPolygon(obj.traceGrids(side),ids(side+1));
          if project
            clip = (clip - origin)*basis;
          end

          % clip the polygons (in 2D local frame)
          verts = clipPolygon(verts,clip);
          if size(verts,1) < 3
            break
          end
        end

        if size(verts,1) < 3
          continue
        end

        [area,center] = computePolygonGeometry(verts);
        if area > 0

          % valid polygon, get 3D coordinates
          if project
            verts = origin + verts*basis.';
            center = origin + center*basis.';
          end
          centers(k,:) = center;
          areas(k) = area;
          counts(k) = numel(verts);
          polygons{k} = reshape(verts.',[],1);
        end
      end

      if ~any(areas > 0)
        error('InterfaceGrid:NoRemainingConnections','No remaining connections in interface grid.');
      end
      verts = ArrayOfArrays(vertcat(polygons{:}),counts);
      obj.intersections = struct('center',centers,'area',areas,'verts',verts);
      intersections = obj.intersections;

    end


    function verts = getIntersectionPolygon(obj,tId)
      % tId is the row of the triplet in neighbors.
      intersects = obj.findGridIntersections();
      coords = intersects.verts.getRowsMatrix(tId);
      dim = size(intersects.center,2);
      verts = reshape(coords,dim,[]).';
    end

  end


  methods (Static, Access = private)


    function verts = cellPolygon(grid,id)
      nodes = grid.getCellNodes(id);
      if isfield(grid.cells,'VTKType') && grid.cells.VTKType(id) == VTKType.Quad9
        % Grid geometry uses the straight edges between the corner nodes.
        nodes = nodes(1:4);
      end
      verts = grid.coordinates(nodes(nodes > 0),:);
    end

  end
end
