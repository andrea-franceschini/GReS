function DruckerPrager_analytical(mat, z0, gamma, M1, M2, numerical_stress_history, elem_coords)
    % Analytical solution for Drucker-Prager Law
    % 
    % INPUT:
    % - mat: Materials object
    % - topology: mesh topology
    % - z0, gamma, M1, M2: parameters for initial stress state
    % - numerical_stress_history: cell array {n_elements}{n_timesteps}(6x1) with numerical stress history
    % - elem_coords: coordinate z degli elementi [n_elements x 1]
    
    % Number of elements
    nz = length(elem_coords);  
    z = elem_coords;  
    
    % Get Material parameters from materials class
    E = mat.getMaterial(1).ConstLaw.E;          
    nu = mat.getMaterial(1).ConstLaw.nu;
    h = mat.getMaterial(1).ConstLaw.h;
    phi = mat.getMaterial(1).ConstLaw.phi;       
    psi = mat.getMaterial(1).ConstLaw.psi;
    
    phi_rad = deg2rad(phi);
    psi_rad = deg2rad(psi);
    
    % Drucker-Prager parameters
    alpha = (3*tan(phi_rad)) / sqrt(9 + 12*tan(phi_rad)^2)/3;
    beta = (3*tan(psi_rad)) / sqrt(9 + 12*tan(psi_rad)^2)/3;  
    xi = 3./(sqrt(9+12*tan(phi_rad).^2));
    H = (h*xi^2)/((alpha+1/sqrt(3))*sqrt(1/3+2*beta^2));
    
    % Elastic constants
    G = E / (2*(1+nu));                      
    K = E / (3*(1-2*nu));                     

    % Hardening parameter
    j = H * (alpha + 1/sqrt(3)) * sqrt(1/3 + 2*beta^2);

    % Time
    t = 0.5:0.5:10;
    n_steps = length(t);
    
    % Initial stress state
    sigmaz = gamma * z0;
    sigmax = M1 * sigmaz;
    sigmay = M2 * sigmaz;
    sigma_initial = zeros(nz, 6);
    sigma_initial(:, 1) = sigmax;
    sigma_initial(:, 2) = sigmay;
    sigma_initial(:, 3) = sigmaz;
    
    stress_path = zeros(nz, 6, n_steps);
    epsilon_total = zeros(nz, 6, n_steps);
    S_history = zeros(nz, n_steps);
    p_history = zeros(nz, n_steps);
    vol_strain_history = zeros(nz, n_steps);
    epsilon_initial = zeros(nz, 6);
    
    % For each element
    for elem_idx = 1:nz
        fprintf('Processing element %d/%d at z = %.3f\n', elem_idx, nz, z(elem_idx));
        
        elem_stress_history = numerical_stress_history{elem_idx};
       
        % Initial tension
        sigma_n = sigma_initial(elem_idx, :)';
        epsilon_n = epsilon_initial(elem_idx, :)';
        
        % Initial pressure
        p_n = sum(sigma_n(1:3)) / 3;
        
        s_n = sigma_n;
        s_n(1:3) = sigma_n(1:3) - p_n;
        
        % S = ||s||
        S_n = sqrt(s_n(1)^2 + s_n(2)^2 + s_n(3)^2 + 2*(s_n(4)^2 + s_n(5)^2 + s_n(6)^2));
        
        % Solution at each step
        for step = 1:n_steps
            if step == 1
                dt = t(step);
                sigma_prev = sigma_n;
            else
                dt = t(step) - t(step-1);
                sigma_prev = elem_stress_history{step-1};
            end
            
            sigma_current = elem_stress_history{step};
            
            sigma_dot = (sigma_current - sigma_prev) / dt;
            
            p_dot = sum(sigma_dot(1:3)) / 3;
            s_dot = sigma_dot;
            s_dot(1:3) = sigma_dot(1:3) - p_dot;
            S_dot = sqrt(s_dot(1)^2 + s_dot(2)^2 + s_dot(3)^2 + 2*(s_dot(4)^2 + s_dot(5)^2 + s_dot(6)^2));
            
            if S_n > 1e-12 && S_dot > 1e-12
                cos_omega_n = (s_n' * s_dot) / (S_n * S_dot);
                cos_omega_n = max(-1, min(1, cos_omega_n));
                omega_n = acos(cos_omega_n);
                % Check plastic loading
%                 if p_dot > -(S_dot/(3*sqrt(2)*alpha)) * cos_omega_n
%                      fprintf("Plastic Loading - element %d, step %d\n", elem_idx, step);
%                 end

                % Use a specific solver (general or deviatoric radial loading)
                if abs(omega_n) < 1e-5
                    [sol] = solve_radial_step(s_n, S_n, epsilon_n, dt, omega_n, S_dot, ...
                                      alpha, beta, j, p_dot, K, G, s_dot);
                elseif abs(omega_n - pi) < 1e-5
                    [sol] = solve_radial_step(s_n, S_n, epsilon_n, dt, omega_n, S_dot, ...
                                      alpha, beta, j, p_dot, K, G, s_dot);
                else
                    [sol] = solve_general_step(s_n, S_n, epsilon_n, dt, omega_n, S_dot, ...
                                              alpha, beta, j, p_dot, K, G, s_dot, elem_idx);
                end
            end

            epsilon_n = sol.epsilon;
            S_n = sol.S;
            
            p_n = p_n + p_dot * dt;
            s_n = s_n + s_dot * dt;
            
            % Save results
            stress_path(elem_idx, :, step) = sigma_current';
            epsilon_total(elem_idx, :, step) = epsilon_n';
            S_history(elem_idx, step) = S_n;
            p_history(elem_idx, step) = p_n;
            vol_strain_history(elem_idx, step) = epsilon_n(1) + epsilon_n(2) + epsilon_n(3);
        end
    end
    
    % Calculate displacement
    uz = zeros(nz, n_steps);
    for step = 1:n_steps
        if nz > 1
            uz(:, step) = cumtrapz(z, epsilon_total(:, 3, step));
        else
            uz(:, step) = epsilon_total(:, 3, step) .* z;
        end
    end
    
    save('DruckerPrager_Analytical.mat', "uz", "z", "t", "sigma_initial", ...
         "stress_path", "epsilon_total", "S_history", "p_history", "vol_strain_history");
    
    fprintf('Done Computing Drucker Prager Analytical solution.\n');
end

function [result] = solve_radial_step(s_n, S_n, epsilon_n, dt, omega_n, S_dot, alpha, beta, j, p_dot, K, G, s_dot)
    epsilon_vol_n = epsilon_n(1) + epsilon_n(2) + epsilon_n(3);
    
    e_n = epsilon_n;
    e_n(1:3) = epsilon_n(1:3) - epsilon_vol_n/3;
    
    tolerance = 1e-6;
    if abs(S_dot) < tolerance
        q = 0;
    elseif abs(omega_n) < tolerance
        q = 1;
    elseif abs(omega_n - pi) < tolerance
        q = -1;
    end
    
    S = S_n + q * S_dot * dt;
    
    % Deviatoric strain
    term1 = zeros(6, 1);
    term2 = zeros(6, 1);
    
    if S_n > 1e-12
        term1 = ((3 * alpha * p_dot * dt) / (sqrt(2) * S_n * j)) * s_n;
    end
    
    coeff = (1/(2*G) + q/(2*j)) * dt;
    term2 = coeff * s_dot;
    
    e = e_n + term1 + term2;
    
    % Volumetric strain
    volumetric_term1 = (1/(3*K) + 3*alpha*beta/j) * p_dot * dt;
    volumetric_term2 = q * beta * S_dot * dt / (sqrt(2) * j);
    
    epsilon_vol = epsilon_vol_n + volumetric_term1 + volumetric_term2;
    
    % Total strain
    epsilon = e;
    epsilon(1:3) = e(1:3) + epsilon_vol/3;
    
    result.epsilon = epsilon;
    result.deviatoric_strain = e;
    result.volumetric_strain = epsilon_vol;
    result.S = S;
end

function [result] = solve_general_step(s_n, S_n, epsilon_n, dt, omega_n, S_dot_norm, alpha, beta, j, p_dot, K, G, s_dot, elem_idx)
    % Solution for general loading
    
    epsilon_vol_n = epsilon_n(1) + epsilon_n(2) + epsilon_n(3);
    
    e_n = epsilon_n;
    e_n(1:3) = epsilon_n(1:3) - epsilon_vol_n/3;
    
    denominator = S_n * cos(omega_n) + S_dot_norm * dt;
    if abs(denominator) < 1e-12
        omega = pi/2;  
    else
        omega = atan((S_n * sin(omega_n)) / denominator);
        if omega < 0
            omega = omega + pi;
        end
    end
    
    if abs(sin(omega)) < 1e-12
        S = S_n;  
    else
        S = S_n * sin(omega_n) / sin(omega);
    end
    
    tr_sigma_dot = sum(s_dot(1:3)) + 3*p_dot;  
    
    if abs(S) > 1e-12 && abs(S_n) > 1e-12
        A = (1/(2*j)) * log(S/S_n);
        if abs(S_dot_norm) > 1e-12
            A = A + (alpha * tr_sigma_dot) / (sqrt(2) * j * S_dot_norm) * log(abs(tan(omega_n/2) / tan(omega/2)));
        end
    else
        A = 0;
    end
    
    B_term1 = dt * (1/G + 1/j) / 2;
    
    B_term2 = 0;
    if abs(S_dot_norm) > 1e-12
        B_term2 = (alpha * tr_sigma_dot * (S - S_n)) / (sqrt(2) * j * S_dot_norm^2);
    end
    
    B_term3 = 0;
    if abs(S_dot_norm) > 1e-12 && abs(tan(omega_n)) > 1e-12
        B_term3 = -(S_n * sin(omega_n) / S_dot_norm) * (A / tan(omega_n) + (omega_n - omega) / (2*j));
    end
    
    B = B_term1 + B_term2 + B_term3;
    
    % Deviatoric strain
    e = e_n + A * s_n + B * s_dot;
    
    % Volumetric strain
    volumetric_term1 = (1/(3*K) + 3*alpha*beta/j) * p_dot * dt;
    
    volumetric_term2 = beta * (S - S_n) / (sqrt(2) * j);
    
    epsilon_vol = epsilon_vol_n + volumetric_term1 + volumetric_term2;
    
    % Total strain
    epsilon = e;
    epsilon(1:3) = e(1:3) + epsilon_vol/3;
    result.epsilon = epsilon;
    result.deviatoric_strain = e;
    result.volumetric_strain = epsilon_vol;
    result.S = S;
    result.omega_final = omega;
end