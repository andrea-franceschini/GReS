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

%% PLOT RESULTS OF MECHANICAL PATCH TEST 2
ms = 3.5;
lwThick = 2;
lwThin = 0.35;
out = load('out1.mat');
figure(1)
hold on
x = out.conforming.x;
tx = out.conforming.tx;
[x,tx] = reOrderTraction(x,tx);
plot(x, tx, 'k-', 'LineWidth', lwThick)
x = out.P0.x;
tx = out.P0.tx;
plot(x, tx, 'g-o', 'LineWidth', lwThin,'MarkerSize',ms,'MarkerEdgeColor','g','MarkerFaceColor','g')
x = out.dual.x;
tx = out.dual.tx;
[x,tx] = reOrderTraction(x,tx);
plot(x, tx, 'r-o', 'LineWidth', lwThin,'MarkerSize',ms,'MarkerEdgeColor','r','MarkerFaceColor','r')
x = out.standard.x;
tx = out.standard.tx;
[x,tx] = reOrderTraction(x,tx);
plot(x, tx, 'b-o', 'LineWidth', lwThin,'MarkerSize',ms,'MarkerEdgeColor','b','MarkerFaceColor','b')

%legend('Standard Conforming', 'P0', 'Dual', 'Standard', 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$t_x$', 'Interpreter', 'latex', 'FontSize', 12);

% Set axis numbers to LaTeX font and increase font size
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 15;

exportgraphics(figure(1), 'Figures/patchTest2_tx1.pdf', 'Resolution', 300);

% ty

figure(2)
hold on
x = out.conforming.x;
ty = out.conforming.ty;
[x,ty] = reOrderTraction(x,ty);
plot(x, ty, 'k-', 'LineWidth', lwThick)
x = out.P0.x;
ty = out.P0.ty;
plot(x, ty, 'g-o', 'LineWidth', lwThin,'MarkerSize',ms,'MarkerEdgeColor','g','MarkerFaceColor','g')
x = out.dual.x;
ty = out.dual.ty;
[x,ty] = reOrderTraction(x,ty);
plot(x, ty, 'r-o', 'LineWidth', lwThin,'MarkerSize',ms,'MarkerEdgeColor','r','MarkerFaceColor','r')
x = out.standard.x;
ty = out.standard.ty;
[x,ty] = reOrderTraction(x,ty);
plot(x, ty, 'b-o', 'LineWidth', lwThin,'MarkerSize',ms,'MarkerEdgeColor','b','MarkerFaceColor','b')

%legend('P0', 'Standard Conforming', 'Dual', 'Standard', 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$t_y$', 'Interpreter', 'latex', 'FontSize', 12);

% Set axis numbers to LaTeX font and increase font size
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 15;

exportgraphics(figure(2), 'Figures/patchTest2_ty1.pdf', 'Resolution', 300);


function [x,t] = reOrderTraction(x,t)
   % reorder the traction array and x-coordinates 
   % to fix bug due to gmsh ordering of cells
   d = diff(x);
   dx = d(1);
   x = x(1):dx:x(end-1);
   tEnd = t(end);
   t(21:39) = t(20:38);
   t(20) = tEnd;
end