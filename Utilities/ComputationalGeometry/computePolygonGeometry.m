function [elementCenter, elementArea] = computePolygonGeometry(poly, normal)

% mex accelerated computational geometry for 3D polygons

idx = mxOrderPointsCCW(poly, normal);
polyCCW = poly(idx,:);

elementCenter = mxComputePolygonCentroid3D(polyCCW);

elementArea = computePolygonArea(polyCCW);

end

