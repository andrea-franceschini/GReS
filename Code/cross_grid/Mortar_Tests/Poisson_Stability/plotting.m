x = load("xCoord.mat");
x_p0 = x.P0;
x_m1 = x.m1;
x_m2 = x.m2;
m1 = load("method1.mat");
m2 = load("method2.mat");
P0 = m1.P0;
dual_1 = m1.DUAL;
standard_1 = m1.STANDARD;
dual_2 = m2.DUAL;
standard_2 = m2.STANDARD;

% First Figure
figure(1)
plot(x_p0, P0, 'k-s', 'LineWidth', 1)
hold on
plot(x_m1, dual_1, 'r-s', 'LineWidth', 1)
plot(x_m1, standard_1, 'b-s', 'LineWidth', 1)

ylim([-2 5]);

legend('Constant', 'Dual', 'Standard', 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$t$', 'Interpreter', 'latex', 'FontSize', 12);

% Set axis numbers to LaTeX font and increase font size
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 12;

% Second Figure
figure(2)
plot(x_m2, dual_2, 'r-s', 'LineWidth', 1)
hold on
plot(x_m2, standard_2, 'b-s', 'LineWidth', 1)

ylim([-2 5]);

legend('Dual', 'Standard', 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$t$', 'Interpreter', 'latex', 'FontSize', 12);

% Set axis numbers to LaTeX font and increase font size
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 12;

%% Third Figure - finer slave - instabilities
figure(3)
xCoord2 = load('xCoord2.mat');
out = load('out.mat');
x_p0 = xCoord2.p0;
x_col = xCoord2.colocated;
dual = out.dual;
standard = out.standard;
p0 = out.p0;
%plot(x_p0, p0, 'k-s', 'LineWidth', 1)
plot(x_col, standard, 'b-s', 'LineWidth', 1)
hold on
plot(x_col, dual, 'r-s', 'LineWidth', 1)

ylim([-3 3]);

legend('Dual', 'Standard', 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$t$', 'Interpreter', 'latex', 'FontSize', 12);

% Set axis numbers to LaTeX font and increase font size
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 12;
exportgraphics(figure(3), 'oscillation_refinment.pdf', 'Resolution', 300);
%% fig. 4 - multipliers convergence profile
figure(4)
out = load('outConv.mat');
h = out.h(1:end-1);
%dual = out.L2dual;
%standard = out.L2stand;
loglog(h, dual, 'k-s', 'LineWidth', 1,'MarkerSize',8)
hold on
loglog(h, standard, 'k-^', 'LineWidth', 1,'MarkerSize',8)

legend('Dual', 'Standard', 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
xlabel('$h_1$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\|t_h - t\|_{-1/2,\Gamma_{f_1}}$', 'Interpreter', 'latex', 'FontSize', 12);

% Set axis numbers to LaTeX font and increase font size
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 12;

% Slightly increase axis limits
xlim([min(h) * 0.9, max(h) * 1.1]);  % Expanding x limits
ylim([min([dual; standard]) * 0.9, max([dual; standard]) * 1.1]);  % Expanding y limits

% Export high-resolution figure
% exportgraphics(figure(1), 'method1.pdf', 'Resolution', 300);
% exportgraphics(figure(2), 'method2.pdf', 'Resolution', 300);
exportgraphics(figure(4), 'conv_profile.pdf', 'Resolution', 300);

%% fig 5: non integer resolution ratio
figure(5)
xCoord2 = load('xCoord3.mat');
out = load('out2.mat');
x_p0 = xCoord2.p0;
x_col = xCoord2.colocated;
dual = out.dual;
standard = out.standard;
p0 = out.p0;
plot(x_p0, p0, 'k-s', 'LineWidth', 1)
hold on
plot(x_col, dual, 'r-s', 'LineWidth', 1)
plot(x_col, standard, 'b-s', 'LineWidth', 1)

ylim([-3 3]);

legend('Constant', 'Dual', 'Standard', 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$t$', 'Interpreter', 'latex', 'FontSize', 12);

% Set axis numbers to LaTeX font and increase font size
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 12;
exportgraphics(figure(5), 'oscillation_nonInteger.pdf', 'Resolution', 300);