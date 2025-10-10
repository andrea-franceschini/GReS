% compare_polygon_intersections.m
% Compares polyshape (default), polyshape (no simplify), and polybool

clc; clear;

% Number of random quadrilateral pairs
N = 1000;

% Domain limits for points
xlim = [-10, 10];
ylim = [-10, 10];

% Preallocate random polygons
X1 = rand(N, 4) * diff(xlim) + xlim(1);
Y1 = rand(N, 4) * diff(ylim) + ylim(1);
X2 = rand(N, 4) * diff(xlim) + xlim(1);
Y2 = rand(N, 4) * diff(ylim) + ylim(1);

% --- Method 1: polyshape (default) ---
t1 = tic;
for i = 1:N
    K1 = convhull(X1(i,:), Y1(i,:));
    K2 = convhull(X2(i,:), Y2(i,:));
    pg1 = polyshape(X1(i,K1), Y1(i,K1));
    pg2 = polyshape(X2(i,K2), Y2(i,K2));
    inter = intersect(pg1, pg2); %#ok<NASGU>
end
time_polyshape_default = toc(t1);

% --- Method 2: polyshape with 'Simplify', false ---
t2 = tic;
for i = 1:N
    K1 = convhull(X1(i,:), Y1(i,:));
    K2 = convhull(X2(i,:), Y2(i,:));
    pg1 = polyshape(X1(i,K1), Y1(i,K1), 'Simplify', false);
    pg2 = polyshape(X2(i,K2), Y2(i,K2), 'Simplify', false);
    inter = intersect(pg1, pg2); %#ok<NASGU>
end
time_polyshape_nosimplify = toc(t2);

% --- Method 3: polycip mex (deprecated) ---
t3 = tic;
for i = 1:N
    K1 = convhull(X1(i,:), Y1(i,:));
    K2 = convhull(X2(i,:), Y2(i,:));
    [x,y] = polyclip(X1(i,K1)', Y1(i,K1)', X2(i,K2)', Y2(i,K2)',1);
end
time_polybool = toc(t3);

% --- Report results ---
fprintf('\n=== Polygon Intersection Timing (N = %d) ===\n', N);
fprintf('polyshape (default)        : %.4f sec\n', time_polyshape_default);
fprintf('polyshape (no simplify)    : %.4f sec\n', time_polyshape_nosimplify);
fprintf('polyclip (mex-file)      : %.4f sec\n', time_polybool);

fprintf('\nAverage time per intersection:\n');
fprintf('polyshape (default)        : %.6f sec\n', time_polyshape_default/N);
fprintf('polyshape (no simplify)    : %.6f sec\n', time_polyshape_nosimplify/N);
fprintf('polybool (mex-file)      : %.6f sec\n', time_polybool/N);


%%

% Random number of vertices (between 3 and 10)
n1 = randi([3, 10]);
n2 = randi([3, 10]);

% Generate random points for each polygon
pts1 = rand(n1, 2) * 10;   % Scale up for visibility
pts2 = rand(n2, 2) * 10;

% Compute convex hull to order vertices (ensure non-intersecting polygon)
K1 = convhull(pts1(:,1), pts1(:,2));
K2 = convhull(pts2(:,1), pts2(:,2));

% Final polygons
c1 = pts1(K1, :);
c2 = pts2(K2, :);

[cx,cy]=polyclip(c1(:,1),c1(:,2),c2(:,1),c2(:,2),1);
p1 = polyshape(c1(:,1),c1(:,2));
p2 = polyshape(c2(:,1),c2(:,2));
pint = polyshape(cx{:},cy{:});
figure(1)
plot(p1, 'FaceColor', 'red', 'FaceAlpha', 0.4, 'LineWidth', 1.5);
hold on
plot(p2,  'FaceColor', 'blue', 'FaceAlpha', 0.4, 'LineWidth', 1.5);
plot(pint,    'FaceColor', 'green', 'FaceAlpha', 0.8, 'LineWidth', 2);
