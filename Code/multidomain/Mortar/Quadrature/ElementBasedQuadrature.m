classdef ElementBasedQuadrature < MortarQuadrature
  % Assign the interface facet's Gauss points to supported trace triplets.

  properties (Access = private)
    remainingPoints
  end

  methods
    function obj = ElementBasedQuadrature(multType,grid,input)

      obj@MortarQuadrature(multType,grid,input);
      
    end

    function [gpCoords,gpWeights] = processMortarPair(obj,triplet,elems)

      ids = obj.grid.neighbors(triplet,:);
      elId = ids(1);
      elem = elems{1};
      available = obj.remainingPoints{elId};
      if isempty(available)
        % initialize list of available gauss points
        available = true(elem.getNumbGaussPts(),1);
      end

      active = find(available);
      gpCoords = {zeros(0,2),zeros(0,2),zeros(0,2)};
      gpWeights = zeros(0,1);
      if isempty(active)
        % no gauss points remained to project
        return
      end

      rule = elem.getGauss();
      coords = elem.getElementCoords(elId);
      points = elem.getGPointsLocation(elId);
      origin = coords(1,:);
      basis = obj.facetBasis(coords);
      points = points(active,:) - origin;
      normals = [];
      if size(coords,2) == 3 && isfield(obj.grid.cells,'avgNodNormal')
        nodes = obj.grid.getCellNodes(elId);
        nodalNormals = obj.grid.cells.avgNodNormal(nodes(1:elem.nNode),:);
        normals = elem.computeBasisF(rule.coord(active,:))*nodalNormals;
      end

      [xiSlave,inSlave] = obj.projectGP(elems{2},ids(2),points,origin,basis,normals);
      [xiMaster,inMaster] = obj.projectGP(elems{3},ids(3),points,origin,basis,normals);
      
      keep = inSlave & inMaster;
      chosen = active(keep);
      gpCoords = {rule.coord(chosen,:),xiSlave(keep,:),xiMaster(keep,:)};
      weights = elem.getDerBasisFAndDet(elId);
      gpWeights = sign(sum(weights))*reshape(weights(chosen),[],1);

      % A point on a shared trace edge can pass both support tests.
      % Consume it on the first matching triplet to integrate it only once.
      available(chosen) = false;
      obj.remainingPoints{elId} = available;
    end
  end

  methods (Access = protected)
    function prepareQuadrature(obj)
      obj.remainingPoints = cell(obj.grid.cells.num,1);
    end
  end

  methods (Access = private)

    function [xi,inside] = projectGP(obj,elem,id,points,origin,basis,normals)
      coords = elem.getElementCoords(id) - origin;
      if isempty(normals)
        [xi,inside] = obj.mapToReference(elem,coords*basis,points*basis);
      else
        xi = zeros(size(points,1),2);
        inside = false(size(points,1),1);
        for g = 1:size(points,1)
          frame = basis;
          normalLength = norm(normals(g,:));
          if isfinite(normalLength) && normalLength > eps
            % Dropping the normal coordinate gives the projection along
            % the interpolated nodal normal before the Newton solve.
            [~,d1,d2] = pointToSurfaceProjection(zeros(1,3), ...
              normals(g,:)/normalLength,points(g,:));
            frame = [d1,d2];
          end
          [xi(g,:),inside(g)] = obj.mapToReference(elem,coords*frame,points(g,:)*frame);
        end
      end
      inside = inside & elem.checkInRange(xi);
    end
  end
end
