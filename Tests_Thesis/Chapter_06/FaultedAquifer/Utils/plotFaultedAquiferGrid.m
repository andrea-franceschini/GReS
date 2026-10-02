function plotFaultedAquiferGrid(grid, info, config)
% Geometry and nonmatching fault sides for comparison with Figure 6.10.
fig = figure('Color', 'w', 'Position', [100 100 1200 520]);
tiledlayout(1, 2, 'TileSpacing', 'compact');
colors = [0.60 0.43 0.35; 0.74 0.62 0.55; 0.68 0.70 0.28];
nexttile;
hold on;
faces = find(grid.surfaces.tag >= 3);
owner = grid.faces.neighbors(grid.surfaces.faceId(faces), 1);
patch('Vertices', grid.coordinates, 'Faces', grid.surfaces.connectivity(faces, :), ...
    'FaceVertexCData', grid.cells.tag(owner), 'FaceColor', 'flat', ...
    'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.15);
for k = 1:2
    xyz = grid.cells.center(info.wellCells(k), :);
    top = max(grid.coordinates(:, 3)) + 3;
    plot3([xyz(1) xyz(1)], [xyz(2) xyz(2)], [xyz(3) top], 'r-', 'LineWidth', 1.5);
    plot3(xyz(1), xyz(2), top, 'ro', 'MarkerFaceColor', 'r');
end
axis equal tight;
view(35, 24);
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('Native faulted aquifer mesh');
colormap(colors);
clim([0.5 3.5]);
bar = colorbar;
bar.Ticks = 1:3;
bar.TickLabels = info.materialNames;
nexttile;
hold on;
for tag = [1 2]
    faces = find(grid.surfaces.tag == tag);
    patch('Vertices', grid.coordinates, 'Faces', grid.surfaces.connectivity(faces, :), ...
        'FaceColor', 'none', 'EdgeColor', colors(2 * tag - 1, :), 'LineWidth', 0.65);
end
axis equal tight;
view(90, 0);
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title(sprintf('Nonmatching fault meshes: inclination %g degrees', config.faultAngle));
exportgraphics(fig, 'Outputs/mesh.png', 'Resolution', 200);
exportgraphics(fig, 'Outputs/mesh.pdf', 'ContentType', 'vector');
end
