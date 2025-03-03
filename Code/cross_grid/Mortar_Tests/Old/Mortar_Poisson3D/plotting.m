clear
close all

%% 3D POISSON PROBLEM
ms = 9;
lw = 1.5;
L2eb = load("Results\slaveL2_eb.dat");
L2g4 = load("Results\slaveL2_gauss_Int4.dat");
%L2g6 = load("Results_lin\L2_gauss_Int4.dat");
h = [0.25;0.125;0.0625;0.03125];
loglog(h,L2eb,'k-','LineWidth',lw,'Marker','s','MarkerSize',ms)
hold on
loglog(h,L2g4,'k-','LineWidth',lw,'Marker','^','MarkerSize',ms+1)
%loglog(h,L2g6,'k-','LineWidth',lw,'Marker','^','MarkerSize',ms)
xlabel('h')
ylabel('Error norm')
xlim([2.5e-2 3e-1])
ylim([4e-4 5e-1])
% set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 12);
% ax = get(gca,'XTickLabel');
% set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 12)
%axis tight

H1eb = load("Results\slaveH1_eb.dat");
H1g4 = load("Results\slaveH1_gauss_Int4.dat");
loglog(h,H1eb,'k-','LineWidth',lw,'Marker','s','MarkerSize',ms)
hold on
loglog(h,H1g4,'k-','LineWidth',lw,'Marker','^','MarkerSize',ms+1)

% Set the font name and size for all text in the tiled layout
% allAxes = findall(t, 'Type', 'Axes');
% for i = 1:length(allAxes)
%     set(allAxes(i), 'FontName', 'Liberation Serif', 'FontSize', 12); % Change 'Arial' to your desired font and 14 to your desired size
%     % Update title, xlabel, and ylabel specifically
%     % set(get(allAxes(i), 'XLabel'), 'FontName', 'Times', 'FontSize', 12);
%     % set(get(allAxes(i), 'YLabel'), 'FontName', 'Times', 'FontSize', 12);
% end
legend('Element based', 'RBF Gaussian - M = 4');
%leg.Layout.Tile = 'north';
legend boxoff    
legend('Orientation','vertical')
legend('Location','northwest')
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 15);
ax = get(gca,'XTickLabel');
set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 15)
nameOut = 'Poisson_3D_lin';
stmp = strcat('Plots/', nameOut, '.pdf');
% exportgraphics(gcf,stmp,'Resolution',400)
exportgraphics(gcf,stmp,'ContentType','vector')