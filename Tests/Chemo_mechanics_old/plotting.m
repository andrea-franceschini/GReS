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

% Post-processing - Plotting

%% Figure 1 - Validate concentration profile with Fig.3 from [Zhang,2007]
figure();
hold on;
tol = 1e-1;
plot_idx = find(cx >= -tol & abs(cy) < tol & abs(cz) < tol); % Plotting along x-axis
t_target = 1000; % target time in seconds
t_target_d = t_target*params.D / params.Rp^2;
[~, timestep_idx] = min(abs(output_times - t_target_d));
[cx_sorted, sort_idx] = sort(cx(plot_idx)); % sort x-coordinates
p_sorted = p(timestep_idx, plot_idx(sort_idx));
plot(cx_sorted, params.c_max*p_sorted, '--', 'LineWidth', 3, ...
    'DisplayName', 'GReS model');

load("Validation plots\Zhang2007_5um_1000s_i2A_c.csv");
x_plt = Zhang2007_5um_1000s_i2A_c(:, 1) / params.Rp;
c_plt = Zhang2007_5um_1000s_i2A_c(:, 2);
plot(x_plt, c_plt, 'LineWidth', 3, 'DisplayName', 'COMSOL model');

xlabel('x-coordinate along the x-axis');
ylabel('Li concentration c');
legend('show', 'Location', 'best');
xlim([0 1]);
% exportgraphics(gcf, 'Validation plots/c_Zhang2007_validation.eps', 'Resolution', 300);
% exportgraphics(gcf, 'Validation plots/c_Zhang2007_validation.pdf', 'Resolution', 300);
% exportgraphics(gcf, 'Validation plots/c_Zhang2007_validation.png', 'Resolution', 300);

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

load("Validation plots\Zhang2007_5um_1000s_i2A_u.csv");
x_plt = Zhang2007_5um_1000s_i2A_u(:, 1) / params.Rp;
u_plt = Zhang2007_5um_1000s_i2A_u(:, 2);
plot(x_plt, u_plt, 'LineWidth', 3, 'DisplayName', 'COMSOL model');

xlabel('x-coordinate along the x-axis');
% ylabel('Nondimensional displacement magnitude, |u|');
ylabel('Nondimensional radial displacement, u');
legend('show', 'Location', 'best');
xlim([0 1]);
% ylim([0 0.1]);
% exportgraphics(gcf, 'Validation plots/u_Zhang2007_validation.eps', 'Resolution', 300);
% exportgraphics(gcf, 'Validation plots/u_Zhang2007_validation.pdf', 'Resolution', 300);
% exportgraphics(gcf, 'Validation plots/u_Zhang2007_validation.png', 'Resolution', 300);

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
xlim([0 1]);
% exportgraphics(gcf, 'Validation plots/GReS_sigmar.png', 'Resolution', 300);

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
xlim([0 1]);
% exportgraphics(gcf, 'Validation plots/GReS_sigmat.png', 'Resolution', 300);