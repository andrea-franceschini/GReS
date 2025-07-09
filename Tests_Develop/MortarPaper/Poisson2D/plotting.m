% % Plot parameters
clear
close all
lw = 0.75;  % Line width
ms = 10;    % Marker size
% 
% % Plot
% figure(1); clf
% loglog(h, sb, 'k-', 'LineWidth', 1.5, 'Marker', 's', 'MarkerSize', ms)
% hold on
% loglog(h, eb, 'r-', 'LineWidth', lw, 'Marker', '^', 'MarkerSize', ms);
% loglog(h, rbf_gauss, 'b-', 'LineWidth', lw, 'Marker', '+', 'MarkerSize', ms)
% loglog(h, rbf_wendland, 'b-', 'LineWidth', lw, 'Marker', 'o', 'MarkerSize', ms)
% 
% %legend({'SegmentBased', 'SegmentBased', 'RBF-Gauss $M=6$', 'RBF-Wendland $M=6$'}, 'Location', 'east','Interpreter', 'latex','FontSize',11)
% xlabel('$h_2$', 'Interpreter', 'latex')
% ylabel('$L^2$ error norm', 'Interpreter', 'latex')


%% HEXA

sb = load(strcat('SB','_hexa.mat'));
eb = load(strcat('EB','_hexa.mat'));
rbf = load(strcat('RBF','_hexa_gauss.mat'));

h = 1./(1.5*(2*2.^(1:6-1)));

figure(1); clf
loglog(h, sb.L2, 'k-', 'LineWidth', 1.5, 'Marker', 's', 'MarkerSize', ms)
hold on
loglog(h, eb.L2, 'r-', 'LineWidth', lw, 'Marker', '^', 'MarkerSize', ms);
loglog(h, rbf.L2, 'b-', 'LineWidth', lw, 'Marker', '+', 'MarkerSize', ms)
loglog(h, sb.H1, 'k-', 'LineWidth', 1.5, 'Marker', 's', 'MarkerSize', ms)
loglog(h, eb.H1, 'r-', 'LineWidth', lw, 'Marker', '^', 'MarkerSize', ms);
loglog(h, rbf.H1, 'b-', 'LineWidth', lw, 'Marker', '+', 'MarkerSize', ms)




lims = [0.75*min(h) 1.25*max(h)];
xlim(lims)
ylim([1e-6 1e0]);

xlabel('$h_2$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('Error broken norm', 'Interpreter', 'latex','FontSize', 16)

set(gca, 'TickLabelInterpreter', 'latex','FontSize', 12)
exportgraphics(figure(1),'convLin_2D.pdf','ContentType', 'vector')



%legend({'SegmentBased', 'SegmentBased', 'RBF-Gauss $M=6$', 'RBF-Wendland $M=6$'}, 'Location', 'east','Interpreter', 'latex','FontSize',11)


%% HEXA 27

sb = load(strcat('SB','_hexa27.mat'));
eb = load(strcat('EB','_hexa27.mat'));
rbf = load(strcat('RBF','_hexa27_gauss.mat'));

figure(2); clf
loglog(h, sb.L2, 'k-', 'LineWidth', 1.5, 'Marker', 's', 'MarkerSize', ms)
hold on
loglog(h, eb.L2, 'r-', 'LineWidth', lw, 'Marker', '^', 'MarkerSize', ms);
loglog(h, rbf.L2, 'b-', 'LineWidth', lw, 'Marker', '+', 'MarkerSize', ms)
loglog(h, sb.H1, 'k-', 'LineWidth', 1.5, 'Marker', 's', 'MarkerSize', ms)
loglog(h, eb.H1, 'r-', 'LineWidth', lw, 'Marker', '^', 'MarkerSize', ms);
loglog(h, rbf.H1, 'b-', 'LineWidth', lw, 'Marker', '+', 'MarkerSize', ms)

xlim(lims)
ylim([1e-6 1e0]);

xlabel('$h_2$', 'Interpreter', 'latex','FontSize', 16)
ylabel('Error broken norm', 'Interpreter', 'latex', 'FontSize', 16)

set(gca, 'TickLabelInterpreter', 'latex','FontSize', 12)
exportgraphics(figure(2),'convQuad_2D.pdf','ContentType', 'vector')


