function initializeFaultStress(domain, config)
% Balanced geostatic stress is an initial state, not an extra gravity load.
% Contact.initialize projects it onto the fault and stores the initial traction.
grid = domain.grid;
element = FiniteElementType.create(VTKType.Hexa, grid, 'gaussOrder', 2);
N = element.getBasisFinGPoints();
numGP = element.getNumbGaussPts();
sigma = zeros(numGP * grid.cells.num, 6);
topFaces = find(grid.surfaces.tag == 8);
for cellId = 1:grid.cells.num
    nodes = grid.getCellNodes(cellId);
    points = N * grid.coordinates(nodes, :);
    % Thickness of the vertical native column gives the local overburden.
    center = grid.cells.center(cellId, :);
    [~, nearest] = min(sum((grid.surfaces.center(topFaces, 1:2) - center(1:2)).^2, 2));
    top = grid.surfaces.center(topFaces(nearest), 3);
    depth = max(0, top - points(:, 3));
    vertical = config.specificWeight * depth;
    ids = (cellId - 1) * numGP + (1:numGP);
    sigma(ids, 1) = -config.K0 * vertical;
    sigma(ids, 2) = -config.K0 * vertical;
    sigma(ids, 3) = -vertical;
end
domain.setState(sigma, 'stress');
end
