function [elementCenter, elementArea] = computePolygonGeometry(intersectionPoints, normal)

elementCenter = mean(intersectionPoints,1);

elementArea = computeSurfaceArea(intersectionPoints, normal);

end


function area = computeSurfaceArea(points, normal)

% reorder points CCW
idx = orderPointsCCW(points, normal);
P = points(idx,:);
%P = points;
% fan triangulation
P0 = P(1,:);
area = 0;

for a = 2:size(P,1)-1
    v1 = P(a,:)   - P0;
    v2 = P(a+1,:) - P0;
    area = area + norm(cross(v1, v2));
end

area = 0.5 * area;
end


function indices = orderPointsCCW(points, normal)
% points : N x 3 (coplanar, unordered)
% normal : 1 x 3 (orientation reference)
% indices: permutation such that points(indices,:) are CCW

numPoints= size(points,1);
indices = (1:numPoints).';

% centroid
centroid = mean(points,1);

v0 = centroid - points(1,:);
v0 = v0 / norm(v0);

angle = zeros(numPoints,1);
angle(1) = 0;
for a = 2:numPoints
    v = centroid - points(a,:);
    
    dotv = dot(v, v0);
    
    crossProduct = cross(v, v0);
    detv = dot(normal, crossProduct);
    
    angle(a) = atan2(detv, dotv);
end


[~, perm] = sort(angle);
indices = indices(perm);
end

