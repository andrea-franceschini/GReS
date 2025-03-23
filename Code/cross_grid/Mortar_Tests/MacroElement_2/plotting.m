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

% Create figure
figure(1)
hold on
grid on

% Plot coarse mesh data
plot(log(h), log(p0_coarse), 'r-s', 'LineWidth', 1, 'DisplayName', 'Coarse $p_0$');
plot(log(h), log(dual_coarse), 'k-s', 'LineWidth', 1, 'DisplayName', 'Coarse $\beta_h^*$');

% Plot fine mesh data
plot(log(h), log(p0_fine), 'r--s', 'LineWidth', 1, 'DisplayName', 'Fine $p_0$');
plot(log(h), log(dual_fine), 'k--s', 'LineWidth', 1, 'DisplayName', 'Fine $\beta_h^*$');

% Labels and title with LaTeX
xlabel('$\log(h)$', 'Interpreter', 'latex');
ylabel('$\log(\beta_h^{\star})$', 'Interpreter', 'latex');
%title('Inf-Sup Test Results', 'Interpreter', 'latex');

% Legend with LaTeX
%legend('Interpreter', 'latex', 'Location', 'best');

% Export high-resolution figure
exportgraphics(figure(1), 'infSup_test.pdf', 'Resolution', 300);




