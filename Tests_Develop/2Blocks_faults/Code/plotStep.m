function plotStep(str,tStep)
% Plot solution for a desired time step
% Figure 1: converged stresses and gap along the vertical axis
% Figure 2: convergence profiele of active set loop

% fig = figure(1);
% %title(strcat('Converge profile - Load step_',num2str(tStep)))
% str = outStruct(tStep);
% numbASIter = max(str.itAS);
% k = 0;
% for i = 1:numbASIter
%    id = str.itAS == i;
%    semilogy(k+1:k+sum(id),str.rhsNorm(id),'k-s','LineWidth',1)
%    hold on
%    k = k+sum(id);
% end
% xlabel('Non-Linear Iteration', 'Interpreter', 'latex')
% ylabel('Residual norm', 'Interpreter', 'latex')
% 
% set(gca, 'FontName', 'Times', 'FontSize', 12, 'TickLabelInterpreter', 'latex');
% 
% % Specify figure size (Width x Height in inches)
% width = 7;  % Width in inches
% height = 4; % Height in inches
% set(fig, 'Units', 'Inches', 'Position', [1, 1, width, height]);
% 
% print(strcat('Results/conv_',num2str(tStep),'.png'), '-dpng', '-r500');
% 
% coordZ = linspace(0,15,numel(str.gap));
% min_sn = min((vertcat(outStruct.s_n)));

str = str(tStep);

ns = numel(str.sn);
dz = 15/ns;
z = linspace(dz,15-dz,ns);

fig = figure(1);
t = tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
nexttile
plot(str.sn,z,'k-s','LineWidth',1)
xlabel('$\sigma_n (kPa)$', 'Interpreter', 'latex')
ylabel('z (m)', 'Interpreter', 'latex')
xlim([-7 0])
xticks([-6 -4 -2 0])
ylim([0 15])
set(gca, 'FontName', 'Times', 'FontSize', 12, 'TickLabelInterpreter', 'latex');
%
nexttile
plot(str.t_norm,z,'k-s','LineWidth',1)
xlabel('$\tau_{norm}$', 'Interpreter', 'latex')
ylabel('z (m)', 'Interpreter', 'latex')
xlim([0 5])
xticks([0 2.5 5])
ylim([0 15])
set(gca, 'FontName', 'Times', 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% Specify figure size (Width x Height in inches)
width = 7;  % Width in inches
height = 3; % Height in inches
set(fig, 'Units', 'Inches', 'Position', [1, 1, width, height]);
% export plots
print(strcat('Results/solution_',num2str(tStep),'.png'), '-dpng', '-r500');

end

