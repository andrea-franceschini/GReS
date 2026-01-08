% Note that the function is currently implemented only for 1 Gauss point
% (nG = 1). Strains and stresses are computed using Voigt notation.
function [strain, stress] = computeStrainsAndStresses(output_times, p, ...
    u, mesh, params)
    % Extract the number of timesteps
    nTimesteps = length(output_times);

    % Initialize strain and stress matrices
    strain = zeros(nTimesteps, mesh.nCells, 6);
    stress = zeros(nTimesteps, mesh.nCells, 6);

    % Loop over each cell
    for timestep = 1:nTimesteps
        % Initialize strain_t and stress_t for the current timestep
        t = output_times(timestep);
        strain_t = zeros(mesh.nCells, 6);
        stress_t = zeros(mesh.nCells, 6);
        for el = 1:mesh.nCells
            % Extract the nodes and their coordinates in the cell
            el_nodes = mesh.cells(el, 1:mesh.cellNumVerts(el));
            el_coords = mesh.coordinates(el_nodes, :);
    
            % Get gradN from the cell coordinates
            gradN = computeShapeFunctionGradients(el_coords);
    
            % Compute matrix B (Voigt notation: [exx, eyy, ezz, 2*exy, 2*eyz, 2*exz])
            nNodes = mesh.cellNumVerts(el);
            B = zeros(6, mesh.nDim * nNodes);
            for i = 1:nNodes
                Bi = [ gradN(i,1)  0           0;
                       0            gradN(i,2)  0;
                       0            0           gradN(i,3);
                       gradN(i,2)   gradN(i,1)  0;
                       0            gradN(i,3)  gradN(i,2);
                       gradN(i,3)   0           gradN(i,1)];
                B(:, (mesh.nDim*(i-1)+1):(mesh.nDim*i)) = Bi;
            end
    
            % Extract displacement values for the selected cell
            % Extract DOF indices for mesh.nDim dimensions
            dof = zeros(1, mesh.nDim * mesh.cellNumVerts(el));
            for i = 1:mesh.cellNumVerts(el)
                base = (el_nodes(i)-1)*mesh.nDim;
                dof(mesh.nDim*(i-1)+1:mesh.nDim*i) = base+1:base+mesh.nDim;
            end
            % dof = repelem(mesh.nDim*el_nodes, mesh.nDim); % (4 x nDim)
            % dof = dof + repmat(-(mesh.nDim-1):0,1,mesh.cellNumVerts(el));
    
            % b) Extracting relevant values of displacements
            el_u = u(timestep, dof)';
    
            % Compute strains using strain = B*u
            strain_t(el,:) = B*el_u;
    
            % Extract pressure values for the selected cell
            el_p = p(timestep, el_nodes);

            % Compute stresses based off of strains
            tr_straint_el = sum(strain_t(el,1:3));
            el_p_avg = mean(el_p);
            stress_t(el,:) = 2*params.G_d * strain_t(el,:) + ...
                (params.lambda_d * tr_straint_el - ...
                params.beta_d * (el_p_avg - params.c_in_d))*[1 1 1 0 0 0];
        end
        % Store strain_t and stress_t in strain and stress matrices
        strain(timestep, :, :) = strain_t;
        stress(timestep, :, :) = stress_t;
    end
end