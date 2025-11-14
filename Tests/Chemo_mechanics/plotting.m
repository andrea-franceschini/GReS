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
tol = 1e-1;
plot_idx = find(abs(cy) < tol & abs(cz) < tol); % Plotting along x-axis
t_target = 1000; % target time in seconds
t_target_d = 1000*params.D / params.Rp^2;
[~, timestep_idx] = min(abs(output_times - t_target_d));
[cx_sorted, sort_idx] = sort(cx(plot_idx)); % sort x-coordinates
p_sorted = p(timestep_idx, plot_idx(sort_idx));
plot(cx_sorted, params.c_max*p_sorted, '-o');
xlabel('x-coordinate along the x-axis');
ylabel('Li concentration c');

%% Figure 2 - Check displacements for the same timestep along the x-axis
figure();
hold on;
u_abs_sorted = u_abs(timestep_idx, plot_idx(sort_idx));
plot(cx_sorted, u_abs_sorted);
% for timestep = 1:length(output_times)
%     plot(cx_sorted, u_abs_sorted(timestep, :));
% end
xlabel('x-coordinate along the x-axis');
ylabel('Absolute displacement |u|');
% ylim([0 0.1]);

%% Figure 3 - Check stresses for the same timestep along the x-axis
figure();
hold on;
sigma_xx = stress_nodal(:,:,3);
sigma_xx_sorted = sigma_xx(:, plot_idx(sort_idx));
for timestep = 1:length(output_times)
    plot(cx_sorted, sigma_xx_sorted(timestep, :), 'DisplayName', ...
        num2str(timestep));
end
xlabel('x-coordinate along the x-axis');
ylabel('\sigma_{xx}');
% legend('show', 'location', 'best');