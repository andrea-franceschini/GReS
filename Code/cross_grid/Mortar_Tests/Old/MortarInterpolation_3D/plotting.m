clear
close all

%% INTERPOLATION TEST 1D - linear at left, quadratic at right
ms = 8;
lw = 1.2;
t = tiledlayout(1,2);
nexttile
L2eb = load("Results_lin/L2_eb.dat");
L2g4 = load("Results_lin/L2_gauss_Int4.dat");
L2g6 = load("Results_lin/L2_gauss_Int4.dat");
h = load("Results_quad/h.dat");
loglog(h,L2eb,'k-','LineWidth',lw,'Marker','s','MarkerSize',ms)
hold on
loglog(h,L2g4,'k-','LineWidth',lw,'Marker','o','MarkerSize',ms)
loglog(h,L2g6,'k-','LineWidth',lw,'Marker','^','MarkerSize',ms)
xlabel('h')
ylabel('L2 error norm')
xlim([3e-3 5e-1])
ylim([1e-8 1e0])
% xlim([4 24])
xticks([1e-2 1e-1])
% set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 12);
% ax = get(gca,'XTickLabel');
% set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 12)
%axis tight


nexttile
L2eb = load("Results_quad/L2_eb.dat");
L2g4 = load("Results_quad/L2_gauss_Int4.dat");
L2g6 = load("Results_quad/L2_gauss_Int6.dat");
h = load("Results_quad/h.dat");
loglog(h,L2eb,'k-','LineWidth',lw,'Marker','s','MarkerSize',ms)
hold on
loglog(h,L2g4,'k-','LineWidth',lw,'Marker','o','MarkerSize',ms)
loglog(h,L2g6,'k-','LineWidth',lw,'Marker','^','MarkerSize',ms)
xlim([3e-3 5e-1])
ylim([1e-8 1e0])
% xlim([4 24])
xticks([1e-2 1e-1])
xlabel('h')
ylabel('L2 error norm')

% Set the font name and size for all text in the tiled layout
% allAxes = findall(t, 'Type', 'Axes');
% for i = 1:length(allAxes)
%     set(allAxes(i), 'FontName', 'Liberation Serif', 'FontSize', 12); % Change 'Arial' to your desired font and 14 to your desired size
%     % Update title, xlabel, and ylabel specifically
%     % set(get(allAxes(i), 'XLabel'), 'FontName', 'Times', 'FontSize', 12);
%     % set(get(allAxes(i), 'YLabel'), 'FontName', 'Times', 'FontSize', 12);
% end
% legend('Element based', 'RBF Gaussian - M = 16','RBF Gaussian - M = 36');
% %leg.Layout.Tile = 'north';
% legend boxoff    
% legend('Orientation','vertical')
% legend('Location','northwest')
% set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 12);
% ax = get(gca,'XTickLabel');
% set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 12)
% nameOut = 'Interp_2D';
% stmp = strcat('Plots/', nameOut, '.pdf');
% % exportgraphics(gcf,stmp,'Resolution',400)
% exportgraphics(gcf,stmp,'ContentType','vector')

%%%%

% Set the font name and size for all text in the tiled layout
allAxes = findall(t, 'Type', 'Axes');
for i = 1:length(allAxes)
    set(allAxes(i), 'FontName', 'Liberation Serif', 'FontSize', 12); % Change 'Arial' to your desired font and 14 to your desired size
    % Update title, xlabel, and ylabel specifically
    % set(get(allAxes(i), 'XLabel'), 'FontName', 'Times', 'FontSize', 12);
    % set(get(allAxes(i), 'YLabel'), 'FontName', 'Times', 'FontSize', 12);
end
leg = legend('Element based', 'RBF Gaussian - n_M = 4','RBF Gaussian - n_M = 6');
leg.Layout.Tile = 'north';
legend boxoff    
legend('Orientation','horizontal')
nameOut = 'Interp_2D';
stmp = strcat('Plots/', nameOut, '.pdf');
x0=10;
y0=10;
width=800;
height=400;
set(gcf,'position',[x0,y0,width,height])
% exportgraphics(gcf,stmp,'Resolution',400)
exportgraphics(gcf,stmp,'ContentType','vector')