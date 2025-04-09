%% First Figure - Taylor papadoupoulos patch test
out = load('OutPT1.mat');
x = out.P0.x;
tx = out.P0.tx;
ty = out.P0.ty;
ms = 8;
figure(1)
plot(x, ty, 'k-s', 'LineWidth', 1,'MarkerSize',ms)
hold on
plot(x, ty, 'r-^', 'LineWidth', 1,'MarkerSize',ms)

ylim([ -10.01 -9.99]);

legend('P0 unstabilized', 'P0 stabilized', 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$t_y$', 'Interpreter', 'latex', 'FontSize', 12);

% Set axis numbers to LaTeX font and increase font size
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 12;

exportgraphics(figure(1), 'Figures/patchTest1.pdf', 'Resolution', 300);