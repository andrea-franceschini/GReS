function area = computePolygonArea(points)

% points are in ccw order

% fan triangulation
P0 = points(1,:);
area = 0;

for a = 2:size(points,1)-1
    v1 = points(a,:)   - P0;
    v2 = points(a+1,:) - P0;
    area = area + norm(cross(v1, v2));
end

area = 0.5 * area;
end

