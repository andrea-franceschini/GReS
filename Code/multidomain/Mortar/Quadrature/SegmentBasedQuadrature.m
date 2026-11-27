classdef SegmentBasedQuadrature < MortarQuadrature
  % Integrate over the intersection polygon of each interface triplet.

  properties (Access = private)
    triangle
  end

  methods
    function obj = SegmentBasedQuadrature(multType,grid,input)

      obj@MortarQuadrature(multType,grid,input);
      obj.triangle = Triangle('gaussOrder',obj.gaussOrder);

    end

    function [gpCoords,gpWeights] = processMortarPair(obj,triplet,elems)
      polygon = obj.grid.getIntersectionPolygon(triplet);
      if isempty(polygon)
        gpCoords = {zeros(0,2),zeros(0,2),zeros(0,2)};
        gpWeights = zeros(0,1);
        return
      end

      ids = obj.grid.neighbors(triplet,:);
      interfaceCoords = elems{1}.getElementCoords(ids(1));

      % project polygon in 2D space
      origin = interfaceCoords(1,:);
      basis = obj.facetBasis(interfaceCoords);
      polygon = (polygon - origin)*basis;
      polygon = polygon(mxOrderPointsCCW(polygon),:);

      % The clipped polygon is convex, so a fan covers it without overlap.
      % Store each triangle's points and weights in consecutive rows.
      
      ng = obj.triangle.getNumbGaussPts();
      nTri = size(polygon,1)-2;
      points = zeros(nTri*ng,2);
      gpWeights = zeros(nTri*ng,1);
      for t = 1:nTri
        vertices = polygon([1,t+1,t+2],:);
        rows = (t-1)*ng + (1:ng);
        points(rows,:) = obj.triangle.getGPointsLocation(vertices);
        gpWeights(rows) = obj.triangle.getDerBasisFAndDet(vertices);
      end
      keep = gpWeights ~= 0;
      points = points(keep,:);
      gpWeights = gpWeights(keep);

      % Invert all three reference maps in the same clipping plane.
      gpCoords = cell(1,3);
      for side = 1:3
        coords = (elems{side}.getElementCoords(ids(side)) - origin)*basis;
        [gpCoords{side},converged] = obj.mapToReference(elems{side},coords,points);
        assert(all(converged),'MortarQuadrature:ProjectionFailed', ...
          'Could not map quadrature points for triplet %d, grid %d.',triplet,side);
      end
    end
  end

  methods (Access = protected)
    function prepareQuadrature(obj)
      obj.grid.findGridIntersections();
    end
  end
end
