% Plotting the results of the inf-sup test on the macroelement
h = [1/4; 1/8; 1/16; 1/32; 1/64];

% Load data
s_coarse = load('slaveCoarse.mat');
s_fine = load('slaveFine.mat');

% Extract variables
p0_coarse = s_coarse.p0;
dual_coarse = s_coarse.dual;
p0_fine = s_fine.p0;
dual_fine = s_fine.dual;

p0_fine = zeros(length(p0_fine),1);

% Create figure
figure(1)
hold on
grid on

% Plot coarse mesh data
plot(log(h), p0_coarse, 'r-s', 'LineWidth', 1, 'DisplayName', 'Coarse $p_0$');
plot(log(h), dual_coarse, 'k-s', 'LineWidth', 1, 'DisplayName', 'Coarse $\beta_h^*$');

% Plot fine mesh data
plot(log(h), p0_fine, 'r--s', 'LineWidth', 1, 'DisplayName', 'Fine $p_0$');
plot(log(h), dual_fine, 'k--s', 'LineWidth', 1, 'DisplayName', 'Fine $\beta_h^*$');

% Labels and title with LaTeX
xlabel('$\log(h_1)$', 'Interpreter', 'latex');
ylabel('$\beta_h^{\star}$', 'Interpreter', 'latex');
xlim([-4.5 -1])
ylim([-0.1 0.6])
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
%title('Inf-Sup Test Results', 'Interpreter', 'latex');

% Legend with LaTeX
%legend('Interpreter', 'latex', 'Location', 'best');

% Export high-resolution figure
exportgraphics(figure(1), 'infSup_test.pdf', 'Resolution', 300);




