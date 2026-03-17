function area = computePolygonArea2D(points)

% points are in ccw order

P0 = points(1,:);
area = 0;

for a = 2:size(points,1)-1
    v1 = points(a,:)   - P0;
    v2 = points(a+1,:) - P0;
    area = area + abs(v1(1)*v2(2) - v1(2)*v2(1)); % 2D cross product
end

area = 0.5 * area;
end

