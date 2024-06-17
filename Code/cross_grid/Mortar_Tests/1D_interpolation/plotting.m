clear
close all

%% INTERPOLATION TEST 1D - linear at left, quadratic at right
ms = 7.5;
t = tiledlayout(1,2);
%t.TileSpacing = 'compact';
%t.Padding = 'compact';
nexttile
L2eb = load("Results_lin\L2_eb.dat");
%L2g4 = load("Results_lin\L2_gauss_Int4.dat");
L2g6 = load("Results_lin\L2_gauss_Int6.dat");
%L2g8 = load("Results_lin\L2_gauss_Int8.dat");
L2w6 = load("Results_lin\L2_wendland_Int6.dat");
h = load("Results_lin\h.dat");
loglog(h,L2eb,'k-','LineWidth',1,'Marker','s','MarkerSize',ms)
hold on
%loglog(h,L2g4,'k-','LineWidth',1,'Marker','*','MarkerSize',ms)
loglog(h,L2g6,'k-','LineWidth',1,'Marker','^','MarkerSize',ms)
%loglog(h,L2g8,'k-','LineWidth',1,'Marker','s','MarkerSize',ms)
loglog(h,L2w6,'k-','LineWidth',1,'Marker','o','MarkerSize',ms)
%grid on
xlabel('h')
ylabel('L2 error norm')
ylim([1e-9 1e0])
% xlim([4 24])
xticks([1e-3 1e-2 1e-1])
% set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 12);
% ax = get(gca,'XTickLabel');
% set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 12)
%axis tight


nexttile
L2eb = load("Results\L2_eb.dat");
%L2g4 = load("Results\L2_gauss_Int4.dat");
L2g6 = load("Results\L2_gauss_Int6.dat");
%L2g8 = load("Results\L2_gauss_Int8.dat");
L2w6 = load("Results\L2_wendland_Int6.dat");
h = load("Results\h.dat");
loglog(h,L2eb,'k-','LineWidth',1,'Marker','s','MarkerSize',ms)
hold on
%loglog(h,L2g4,'k-','LineWidth',1,'Marker','*','MarkerSize',ms)
loglog(h,L2g6,'k-','LineWidth',1,'Marker','^','MarkerSize',ms)
%loglog(h,L2g8,'k-','LineWidth',1,'Marker','s','MarkerSize',ms)
loglog(h,L2w6,'k-','LineWidth',1,'Marker','o','MarkerSize',ms)
%grid on
xlabel('h')
ylabel('L2 error norm')
ylim([1e-9 1e0])
xticks([1e-3 1e-2 1e-1])
% set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 12);
% ax = get(gca,'XTickLabel');
% set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 12)

% Set the font name and size for all text in the tiled layout
allAxes = findall(t, 'Type', 'Axes');
for i = 1:length(allAxes)
    set(allAxes(i), 'FontName', 'Liberation Serif', 'FontSize', 12); % Change 'Arial' to your desired font and 14 to your desired size
    % Update title, xlabel, and ylabel specifically
    % set(get(allAxes(i), 'XLabel'), 'FontName', 'Times', 'FontSize', 12);
    % set(get(allAxes(i), 'YLabel'), 'FontName', 'Times', 'FontSize', 12);
end
leg = legend('Element based', 'RBF Gaussian - M = 6','RBF Wendland - M = 6');
leg.Layout.Tile = 'north';
legend boxoff    
legend('Orientation','horizontal')
nameOut = 'Interp_1D';
stmp = strcat('Plots/', nameOut, '.pdf');
% exportgraphics(gcf,stmp,'Resolution',400)
exportgraphics(gcf,stmp,'ContentType','vector')

%%
ms = 8;
lw = 1.2;
L2eb = load("Results_lin\L2_eb.dat");
L2g6 = load("Results_lin\L2_gauss_Int6.dat");
L2w6 = load("Results_lin\L2_wendland_Int6.dat");
h = load("Results\h.dat");
loglog(h,L2eb,'k-','LineWidth',lw,'Marker','s','MarkerSize',ms)
hold on
loglog(h,L2g6,'k-','LineWidth',lw,'Marker','o','MarkerSize',ms)
loglog(h,L2w6,'k-','LineWidth',lw,'Marker','^','MarkerSize',ms)
xlabel('h')
ylabel('L2 error norm')
xlim([1e-4 2e-1])
ylim([1e-9 1e0])
% xlim([4 24])
xticks([1e-3 1e-2 1e-1])
% set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 12);
% ax = get(gca,'XTickLabel');
% set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 12)
%axis tight

L2eb = load("Results\L2_eb.dat");
L2g6 = load("Results\L2_gauss_Int4.dat");
L2w6 = load("Results\L2_wendland_Int6.dat");
loglog(h,L2eb,'k-','LineWidth',lw,'Marker','s','MarkerSize',ms)
hold on
loglog(h,L2g6,'k-','LineWidth',lw,'Marker','o','MarkerSize',ms)
loglog(h,L2w6,'k-','LineWidth',lw,'Marker','^','MarkerSize',ms)

% Set the font name and size for all text in the tiled layout
% allAxes = findall(t, 'Type', 'Axes');
% for i = 1:length(allAxes)
%     set(allAxes(i), 'FontName', 'Liberation Serif', 'FontSize', 12); % Change 'Arial' to your desired font and 14 to your desired size
%     % Update title, xlabel, and ylabel specifically
%     % set(get(allAxes(i), 'XLabel'), 'FontName', 'Times', 'FontSize', 12);
%     % set(get(allAxes(i), 'YLabel'), 'FontName', 'Times', 'FontSize', 12);
% end
legend('Element based', 'RBF Gaussian - M = 6','RBF Wendland - M = 6');
%leg.Layout.Tile = 'north';
legend boxoff    
legend('Orientation','vertical')
legend('Location','northwest')
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 12);
ax = get(gca,'XTickLabel');
set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 12)
nameOut = 'Interp_1D';
stmp = strcat('Plots/', nameOut, '.pdf');
% exportgraphics(gcf,stmp,'Resolution',400)
exportgraphics(gcf,stmp,'ContentType','vector')

