function sigma_node = celltonodeStress(sigma_elem, mesh)
% sigma_elem: nTimesteps x nElem x 6 (Voigt) element stresses, nG = 1
% mesh: struct with .cells (nElem x 4), .coordinates (nNodes x 3)
% returns sigma_node: nNodes x 6

nNodes = size(mesh.coordinates,1);
nElem  = mesh.nCells;
nTimesteps = size(sigma_elem,1);

% Initialize output
sigma_node = zeros(nTimesteps, nNodes, 6);

for t = 1 : nTimesteps
    sigma_node_t = zeros(nNodes,6);
    weight_node = zeros(nNodes,1);

    for e = 1:nElem
        nodes = mesh.cells(e,1:mesh.cellNumVerts(e));   % 1x4
        coords = mesh.coordinates(nodes,:);             % 4x3
        V6 = det([ones(4,1), coords]);                  % 6*volume
        vol = abs(V6)/6;
        w = vol;

        sigma_node_t(nodes,:) = sigma_node_t(nodes,:) + ...
            w * reshape(sigma_elem(t,e,:), 1, []);
        weight_node(nodes)   = weight_node(nodes) + w;
    end
    
    % Avoid divide-by-zero
    nz = weight_node > 0;
    sigma_node_t(nz,:) = sigma_node_t(nz,:) ./ weight_node(nz);

    sigma_node(t,:,:) = sigma_node_t;
end
