function ccw_points = orderPointsCCW2D(points)
% points: N x 3 (coplanar)
% indices: permutation such that points(indices,:) are CCW

% Step 1: Compute centroid
centroid = mean(points, 1);

% Step 2: Compute angles relative to centroid
angles = atan2(points(:,2) - centroid(2), points(:,1) - centroid(1));

% Step 3: Sort points by angle
[~, idx] = sort(angles);
ccw_points = points(idx, :);
end
