function Q = pressureLoadMatrix(grid)
% Q*p = integral(B' * [p,p,p,0,0,0]') over all cells, in nodal DOF order.
% Use GReS's shape functions and quadrature; this helper does not make a mesh.
numCells = grid.cells.num;
rows = zeros(24 * numCells, 1);
cols = rows;
values = rows;
element = FiniteElementType.create(VTKType.Hexa, grid, 'gaussOrder', 2);
for cellId = 1:numCells
    nodes = grid.getCellNodes(cellId);
    coords = grid.coordinates(nodes, :);
    [gradient, dJW] = element.getDerBasisFAndDet(coords);
    assert(all(dJW > 0), 'Nonpositive quadrature Jacobian in the pressure map.');
    B = element.getStrainMatrix(gradient);
    q = sum(pagemtimes(permute(B, [2 1 3]), [1; 1; 1; 0; 0; 0]) .* ...
        reshape(dJW, 1, 1, []), 3);
    range = (cellId - 1) * 24 + (1:24);
    rows(range) = reshape((3 * (nodes(:) - 1) + (1:3))', [], 1);
    cols(range) = cellId;
    values(range) = q(:);
end
Q = sparse(rows, cols, values, 3 * grid.nNodes, numCells);
% A constant pressure must have zero total force in each disconnected block.
% This catches a DOF-order or integration error before running mechanics.
force = reshape(Q * ones(numCells, 1), 3, [])';
assert(max(abs(sum(force, 1))) < 1e-8 * max(1, norm(Q, 1)), ...
    'FaultedAquifer:PressureMap', 'Pressure-to-force map is not balanced.');
end
