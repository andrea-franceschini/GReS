function g = thesisHexMesh(n, limits, order)
% Native GReS mesh; optional in-memory elevation to its Hexa27 ordering.
g = structuredMesh(n(1), n(2), n(3), limits(1, :), limits(2, :), limits(3, :));
if order == 1
    return
end
assert(order == 2, 'Only Q1 and Q2 are supported.');
base = g;
fine = structuredMesh(2 * n(1), 2 * n(2), 2 * n(3), limits(1, :), limits(2, :), limits(3, :));
sz = 2 * n + 1;
[I, J, K] = ndgrid(1:2:2 * n(1), 1:2:2 * n(2), 1:2:2 * n(3));
offset = HexahedronQuadratic.coordLoc + 1;
top = zeros(prod(n), 27);
for a = 1:27
    top(:, a) = sub2ind(sz, I(:) + offset(a, 1), J(:) + offset(a, 2), K(:) + offset(a, 3));
end
% Elevate oriented boundary quads: vertices, four edge midpoints, centre.
c = base.coordinates(base.surfaces.connectivity(:), :);
t = (c - limits(:, 1)') ./ (limits(:, 2) - limits(:, 1))' .* (2 * n) + 1;
t = round(t);
t = reshape(t, size(base.surfaces.connectivity, 1), 4, 3);
surf = zeros(size(t, 1), 9);
for a = 1:4
    v = reshape(t(:, a, :), [], 3);
    surf(:, a) = sub2ind(sz, v(:, 1), v(:, 2), v(:, 3));
    b = mod(a, 4) + 1;
    v = reshape((t(:, a, :) + t(:, b, :)) / 2, [], 3);
    surf(:, 4 + a) = sub2ind(sz, v(:, 1), v(:, 2), v(:, 3));
end
v = reshape(mean(t, 2), [], 3);
surf(:, 9) = sub2ind(sz, v(:, 1), v(:, 2), v(:, 3));
g = Grid();
g.nDim = 3;
g.coordinates = fine.coordinates;
g.cells.connectivity = top;
g.cells.VTKType = 29 * ones(size(top, 1), 1);
g.cells.numVerts = 27 * ones(size(top, 1), 1);
g.cells.tag = ones(size(top, 1), 1);
g.surfaces.connectivity = surf;
g.surfaces.VTKType = 28 * ones(size(surf, 1), 1);
g.surfaces.numVerts = 9 * ones(size(surf, 1), 1);
g.surfaces.tag = base.surfaces.tag;
end
