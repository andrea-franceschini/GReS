% Improved version that aligns with Poromechanics.m approach
% Supports both single-point (centroid) and Gauss point integration
% Note: For linear tetrahedra, gradients are constant, so centroid = any point
function [strain, stress] = computeStrainsAndStresses_improved(output_times, p, ...
    u, mesh, params, elements)
    % Extract the number of timesteps
    nTimesteps = length(output_times);

    % Initialize strain and stress matrices
    % If using Gauss points, we'll average to get cell values
    strain = zeros(nTimesteps, mesh.nCells, 6);
    stress = zeros(nTimesteps, mesh.nCells, 6);

    % Check if elements object is provided (for Gauss point integration)
    useGaussPoints = nargin >= 6 && ~isempty(elements);

    % Loop over each cell
    for timestep = 1:nTimesteps
        strain_t = zeros(mesh.nCells, 6);
        stress_t = zeros(mesh.nCells, 6);
        
        for el = 1:mesh.nCells
            el_nodes = mesh.cells(el, 1:mesh.cellNumVerts(el));
            
            % Extract displacement DOF for the element
            dof = zeros(1, mesh.nDim * mesh.cellNumVerts(el));
            for i = 1:mesh.cellNumVerts(el)
                base = (el_nodes(i) - 1) * mesh.nDim;
                dof(mesh.nDim*(i-1)+1:mesh.nDim*i) = base+1:base+mesh.nDim;
            end
            el_u = u(timestep, dof)';
            
            if useGaussPoints
                % Use Poromechanics approach with Gauss points
                vtkId = mesh.cellVTKType(el);
                elem = elements.getElement(vtkId);
                nG = elem.GaussPts.nNode;
                
                % Get shape function derivatives at Gauss points
                gradN = elem.getDerBasisFAndDet(el, 2); % Returns nG x (nNode*nDim)
                
                % Construct B matrix at each Gauss point (6 x nNode*nDim x nG)
                B = zeros(6, elem.nNode * mesh.nDim, nG);
                B(elem.indB(:,2)) = gradN(elem.indB(:,1));
                
                % Compute strain at each Gauss point
                % B is 6 x (nNode*nDim) x nG, el_u is (nNode*nDim) x 1
                strainGP = reshape(pagemtimes(B, el_u), 6, nG)'; % nG x 6
                
                % Average strain over Gauss points (weighted by volume)
                % For constant strain elements, this is just the mean
                strain_t(el, :) = mean(strainGP, 1);
                
            else
                % Original approach: compute at centroid (single point)
                el_coords = mesh.coordinates(el_nodes, :);
                gradN = computeShapeFunctionGradients(el_coords);
                
                % Compute B matrix (6 x 12 for 4-node tetrahedron)
                B = zeros(6, mesh.nDim * mesh.cellNumVerts(el));
                for i = 1:mesh.cellNumVerts(el)
                    Bi = [gradN(i,1)  0           0;
                          0           gradN(i,2)  0;
                          0           0           gradN(i,3);
                          gradN(i,2)  gradN(i,1)  0;
                          0           gradN(i,3)  gradN(i,2);
                          gradN(i,3)  0           gradN(i,1)];
                    B(:, (mesh.nDim*(i-1)+1):(mesh.nDim*i)) = Bi;
                end
                
                % Compute strain
                strain_t(el, :) = (B * el_u)';
            end
            
            % Extract pressure values
            el_p = p(timestep, el_nodes);
            el_p_avg = mean(el_p);
            
            % Compute stresses (same for both approaches)
            tr_strain_el = sum(strain_t(el, 1:3));
            stress_t(el, :) = 2*params.G_d * strain_t(el, :) + ...
                (params.lambda_d * tr_strain_el - ...
                 params.beta_d * (el_p_avg - params.c_in_d)) * [1 1 1 0 0 0];
        end
        
        % Store results
        strain(timestep, :, :) = strain_t;
        stress(timestep, :, :) = stress_t;
    end
end

