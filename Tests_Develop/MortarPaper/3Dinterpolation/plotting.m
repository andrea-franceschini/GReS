% % Plot parameters
clear
close all
lw = 0.6;  % Line width
ms = 9;    % Marker size
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

dir = 'OUT/HEXA';
sb = load(fullfile(dir,"SegmentBased.mat"));
eb = load(fullfile(dir,"ElementBased.mat"));
rbf_gauss = load(fullfile(dir,'RBF_gauss.mat'));
rbf_gauss_4 = load(fullfile(dir,'RBF_gauss_4.mat'));

h = 1./(1.5*(2*2.^(1:8-1)));

figure(1); clf
loglog(h, sb.L2, 'k-', 'LineWidth', 1.5, 'Marker', 's', 'MarkerSize', ms)
hold on
loglog(h, eb.L2, 'r-', 'LineWidth', lw, 'Marker', '^', 'MarkerSize', ms);
loglog(h, rbf_gauss.L2, 'b-', 'LineWidth', lw, 'Marker', '+', 'MarkerSize', ms)
loglog(h, rbf_gauss_4.L2, 'b-', 'LineWidth', lw, 'Marker', 'o', 'MarkerSize', ms)

lims = [0.75*min(h) 1.25*max(h)];
xlim(lims)
ylim([1e-7 1e0]);

xlabel('$h_2$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$L^2$ error norm', 'Interpreter', 'latex','FontSize', 16)

set(gca, 'TickLabelInterpreter', 'latex','FontSize', 12)
exportgraphics(figure(1),'convLin_2D.pdf','ContentType', 'vector')



%legend({'SegmentBased', 'SegmentBased', 'RBF-Gauss $M=6$', 'RBF-Wendland $M=6$'}, 'Location', 'east','Interpreter', 'latex','FontSize',11)


%% HEXA 27

dir = 'OUT/HEXA27';
sb = load(fullfile(dir,"SegmentBased.mat"));
eb = load(fullfile(dir,"ElementBased.mat"));
rbf_gauss = load(fullfile(dir,'RBF_gauss.mat'));
rbf_gauss_4 = load(fullfile(dir,'RBF_gauss_4.mat'));

h = 1./(1.5*(2*2.^(1:7-1)));

figure(2); clf
loglog(h, sb.L2, 'k-', 'LineWidth', 1.5, 'Marker', 's', 'MarkerSize', ms)
hold on
loglog(h, eb.L2, 'r-', 'LineWidth', lw, 'Marker', '^', 'MarkerSize', ms);
loglog(h, rbf_gauss.L2, 'b-', 'LineWidth', lw, 'Marker', '+', 'MarkerSize', ms)
loglog(h, rbf_gauss_4.L2, 'b-', 'LineWidth', lw, 'Marker', 'o', 'MarkerSize', ms)

xlim(lims)
ylim([1e-7 1e0]);

xlabel('$h_2$', 'Interpreter', 'latex','FontSize', 16)
ylabel('$L^2$ error norm', 'Interpreter', 'latex', 'FontSize', 16)

set(gca, 'TickLabelInterpreter', 'latex','FontSize', 12)
exportgraphics(figure(2),'convQiad_2D.pdf','ContentType', 'vector')




