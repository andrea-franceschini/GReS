function points_ccw = orderPointsCCW3D(points, normal)
% points : N x 3 (coplanar, unordered)
% normal : 1 x 3 (orientation reference)
% indices: permutation such that points(indices,:) are CCW

numPoints= size(points,1);

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
points_ccw = points(perm,:);

end


