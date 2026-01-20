% Post-processing - Decoupling ux, uy, uz from u
nTimesteps = length(output_times);
nNodes = topology.nNodes;
ux = zeros(nTimesteps, nNodes);
uy = zeros(nTimesteps, nNodes);
uz = zeros(nTimesteps, nNodes);
for i = 1:nNodes
    ux(:, i) = u(:, 3*(i-1) + 1);
    uy(:, i) = u(:, 3*(i-1) + 2);
    uz(:, i) = u(:, 3*(i-1) + 3);
end
u_abs = (ux.^2 + uy.^2 + uz.^2).^(0.5);

if ~exist('Results', 'dir')
    mkdir('Results');
end

% Post-processing - Plotting

%% Figure 1 - Validate concentration profile with Fig.3 from [Zhang,2007]
figure();
hold on;
tol = 1e-1;
plot_idx = find(abs(cy) < tol & abs(cz) < tol); % Plotting along x-axis
t_target = 1000; % target time in seconds
t_target_d = t_target*params.D / params.Rp^2;
[~, timestep_idx] = min(abs(output_times - t_target_d));
[cx_sorted, sort_idx] = sort(cx(plot_idx)); % sort x-coordinates
p_sorted = p(timestep_idx, plot_idx(sort_idx));
plot(cx_sorted, params.c_max*p_sorted, '--', 'LineWidth', 3, ...
    'DisplayName', 'GReS model');

xlabel('x-coordinate along the x-axis');
ylabel('Li concentration c');
xlim([-2 2]);
% exportgraphics(gcf, 'Results/c_Zhang2007_validation.png', 'Resolution', 300);

%% Figure 2 - Check displacements for the same timestep along the x-axis
figure();
hold on;
ux_sorted = ux(timestep_idx, plot_idx(sort_idx));
u_abs_sorted = u_abs(timestep_idx, plot_idx(sort_idx));
plot(cx_sorted, u_abs_sorted, '--', 'LineWidth', 3, 'DisplayName', ...
    'GReS model');
% for timestep = 1:length(output_times)
%     plot(cx_sorted, u_abs_sorted(timestep, :));
% end

xlabel('x-coordinate along the x-axis');
% ylabel('Nondimensional displacement magnitude, |u|');
ylabel('Nondimensional radial displacement, u');
xlim([-2 2]);
% ylim([0 0.1]);
% exportgraphics(gcf, 'Results/u_Zhang2007_validation.png', 'Resolution', 300);

%% Figures 3 and 4 - Check stresses for the same timestep along the x-axis
figure();
hold on;
sigma_xx = stress_nodal(:,:,1);
sigma_xx_sorted = sigma_xx(:, plot_idx(sort_idx));
for timestep = 1:length(output_times)
    plot(cx_sorted, sigma_xx_sorted(timestep, :), 'DisplayName', ...
        num2str(timestep));
end
xlabel('x-coordinate along the x-axis');
ylabel('\sigma_{xx}');
% legend('show', 'location', 'best');
xlim([-2 2]);
% exportgraphics(gcf, 'Results/GReS_sigmar.png', 'Resolution', 300);

figure();
hold on;
sigma_yy = stress_nodal(:,:,2);
sigma_yy_sorted = sigma_yy(:, plot_idx(sort_idx));
for timestep = 1:length(output_times)
    plot(cx_sorted, sigma_yy_sorted(timestep, :), 'DisplayName', ...
        num2str(timestep));
end
xlabel('x-coordinate along the x-axis');
ylabel('\sigma_{yy}');
% legend('show', 'location', 'best');
xlim([-2 2]);
% exportgraphics(gcf, 'Results/GReS_sigmat.png', 'Resolution', 300);

%% Comparing contact stresses with diffusion induced stresses
[~, zero_idx] = min(abs(cx_sorted));

% Radial stresses
contact_sigma_r = sigma_xx_sorted(:, zero_idx);
% Find maximum along each row
sigma_xx_noZeroIdx = [sigma_xx_sorted(:, 1:zero_idx-1), ...
    sigma_xx_sorted(:, zero_idx+1:end)];
[max_dis_sigma_r, max_dis_sigma_r_idx] = max(abs(sigma_xx_noZeroIdx), ...
    [], 2);
figure();
hold on;
plot(output_times, max_dis_sigma_r, '--', 'LineWidth', 3, ...
    'DisplayName', 'Maximum absolute DIS');
plot(output_times, abs(contact_sigma_r), 'LineWidth', 3, 'DisplayName', ...
    'Absolute contact stress');
xlabel('Dimensionless time');
ylabel('Maximum absolute radial stresses');
legend('show', 'Location', 'best');
% exportgraphics(gcf, 'Results/contactvsdis_sigmar.png', 'Resolution', 300);

% Tangential stresses
contact_sigma_t = sigma_yy_sorted(:, zero_idx);
% Find maximum along each row
sigma_yy_noZeroIdx = [sigma_yy_sorted(:, 1:zero_idx-1), ...
    sigma_yy_sorted(:, zero_idx+1:end)];
[max_dis_sigma_t, max_dis_sigma_t_idx] = max(sigma_yy_noZeroIdx, ...
    [], 2);
[max_compressive_dis_sigma_t, max_compressive_dis_sigma_t_idx] = ...
    min(sigma_yy_noZeroIdx, [], 2);
figure();
hold on;
plot(output_times, max_dis_sigma_t, 'LineWidth', 3, 'DisplayName', ...
    'Maximum tensile DIS');
plot(output_times, abs(sigma_yy_sorted(:,end)), '--', 'LineWidth', 3, ...
    'DisplayName', 'Maximum compressive DIS');
plot(output_times, abs(contact_sigma_t), 'LineWidth', 3, 'DisplayName', ...
    'Compressive stress at the contact point');
xlabel('Dimensionless time');
ylabel('Quantifying compressive contact stresses');
legend('show', 'Location', 'best');
% exportgraphics(gcf, 'Results/contactvsdis_sigmat.png', 'Resolution', 300);
