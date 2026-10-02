function [grid, info] = generateFaultedAquiferGrid(config)
% Two native structured blocks, mapped to a planar inclined fault.
% Native connectivity is retained. Fault nodes are never merged.
blocks = cell(2, 1);
for side = 1:2
    if side == 1
        counts = config.rockCells;
        % Three coarse rock columns; mild grading towards the fault.
        xAxis = linspace(0, 1, counts(1) + 1).^0.8;
    else
        counts = config.soilCells;
        % Fine near the fault, increasing cell width into the reservoir.
        widths = 0.85.^((counts(1) - 1):-1:0);
        xAxis = [0 cumsum(widths) / sum(widths)];
    end
    yAxis = linspace(0, 1, counts(2) + 1);
    zAxis = linspace(0, 1, counts(3) + 1);
    if side == 2
        % Same outer extent, different interior levels: nonmatching fault mesh.
        zAxis = zAxis + config.faultMeshShift * zAxis .* (1 - zAxis);
    end
    block = structuredMesh(xAxis, yAxis, zAxis);
    param = block.coordinates;
    eta = param(:, 1);
    y = config.width * param(:, 2);
    s = param(:, 3);
    faultWarp = 3 * sin(2 * pi * y / config.width);
    zFault = config.thickness * (s - 1) + faultWarp;
    xFault = config.faultPosition + tand(config.faultAngle) * ...
        (zFault + config.thickness / 2);
    if side == 1
        x = eta .* xFault;
        taper = 1 - eta;
    else
        x = xFault + eta .* (config.length - xFault);
        taper = eta;
    end
    horizonWarp = config.layerWarp * sin(1.5 * pi * x / config.length) .* ...
        cos(pi * y / config.width);
    z = zFault + taper .* horizonWarp;
    block.coordinates = [x y z];

    % Surface tags: fault master=1, slave=2; bottom=3, y-min=4,
    % x-min=5, y-max=6, x-max=7, free top=8.
    oldTags = block.surfaces.tag;
    tags = [3 8 4 6 5 7];
    block.surfaces.tag = reshape(tags(oldTags), [], 1);
    if side == 1
        block.surfaces.tag(oldTags == 6) = 1;
        block.cells.tag(:) = 1;        % rock
    else
        block.surfaces.tag(oldTags == 5) = 2;
        % Material assignment follows logical stratigraphy, not global z.
        logicalCenter = mean(reshape(param(block.cells.connectivity', 3), 8, []), 1)';
        block.cells.tag(:) = 2;        % clay
        block.cells.tag(logicalCenter < config.sandFraction) = 3; % sand
    end
    blocks{side} = block;
end

% Assemble the disconnected blocks as one GReS Grid. Node IDs on the fault
% remain separate even where coordinates coincide; TPFA sees a sealed fault.
grid = Grid();
nodeOffset = size(blocks{1}.coordinates, 1);
grid.coordinates = [blocks{1}.coordinates; blocks{2}.coordinates];
for field = {'connectivity', 'VTKType', 'numVerts', 'tag'}
    key = field{1};
    c = blocks{2}.cells.(key);
    f = blocks{2}.surfaces.(key);
    if strcmp(key, 'connectivity')
        c = c + nodeOffset;
        f = f + nodeOffset;
    end
    grid.cells.(key) = [blocks{1}.cells.(key); c];
    grid.surfaces.(key) = [blocks{1}.surfaces.(key); f];
end
grid.processGeometry();
assert(all(isfinite(grid.cells.volume) & grid.cells.volume > 0), ...
    'FaultedAquifer:InvalidCells', 'The mapped grid contains invalid cells.');

info.rockCells = (1:size(blocks{1}.cells.connectivity, 1))';
info.soilCells = (size(blocks{1}.cells.connectivity, 1) + 1:grid.cells.num)';
info.materialNames = {'Rock', 'Clay', 'Sand'};
info.wellCells = zeros(2, 1);
info.wellTargets = zeros(2, 3);
sand = find(grid.cells.tag == 3);
for k = 1:2
    y = config.wellY(k);
    z = config.thickness * (config.wellDepthFraction - 1) + 3 * sin(2 * pi * y / config.width);
    x = config.faultPosition + tand(config.faultAngle) * (z + config.thickness / 2) + config.wellDistance;
    info.wellTargets(k, :) = [x y z];
    [~, nearest] = min(sum((grid.cells.center(sand, :) - [x y z]).^2, 2));
    info.wellCells(k) = sand(nearest);
end
assert(numel(unique(info.wellCells)) == 2, 'Two distinct sandy well cells are required.');
% Check both fault sides lie in the exact 40-degree plane.
for tag = [1 2]
    ids = unique(grid.surfaces.connectivity(grid.surfaces.tag == tag, :));
    xyz = grid.coordinates(ids, :);
    residual = xyz(:, 1) - config.faultPosition - tand(config.faultAngle) * ...
        (xyz(:, 3) + config.thickness / 2);
    assert(max(abs(residual)) < 1e-8, 'Fault geometry is not planar.');
end
fprintf('Native grid: %d cells, %d nodes; wells: %d, %d.\n', ...
    grid.cells.num, grid.nNodes, info.wellCells);
end
