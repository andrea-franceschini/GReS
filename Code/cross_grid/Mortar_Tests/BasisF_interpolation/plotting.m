%% fig INFLUENCE OF RBF type
N = 4:2:24;
L2_c1 = load("c1_L2");
L2_c2 = load("c2_L2.dat");
nc_c1 = load("c1_cond.dat");
nc_c2 = load("c2_cond.dat");
% gauss imq wendland
t = tiledlayout(2,2);
%t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile
semilogy(N,L2_c1(:,1),'k','LineWidth',1,'Marker','s','MarkerSize',6)
hold on
semilogy(N,L2_c1(:,2),'k','LineWidth',1,'Marker','^','MarkerSize',6)
semilogy(N,L2_c1(:,3),'k','LineWidth',1,'Marker','*','MarkerSize',6)
grid on
xlabel('M')
ylabel('RMSE')
ylim([1e-10 1e-1])
xlim([4 24])
xticks(4:4:24)
yticks([1e-10 1e-5 1e-1])


nexttile
semilogy(N,nc_c1(:,1),'k','LineWidth',1,'Marker','s','MarkerSize',6)
hold on
semilogy(N,nc_c1(:,2),'k','LineWidth',1,'Marker','^','MarkerSize',6)
semilogy(N,nc_c1(:,3),'k','LineWidth',1,'Marker','*','MarkerSize',6)
grid on
xlabel('M')
ylabel('Condition number')
ylim([1e0 1e22])
xlim([4 24])
xticks(4:4:24)
yticks([1e0 1e5 1e10 1e15 1e20])

nexttile
semilogy(N,L2_c2(:,1),'k','LineWidth',1,'Marker','s','MarkerSize',6)
hold on
semilogy(N,L2_c2(:,2),'k','LineWidth',1,'Marker','^','MarkerSize',6)
semilogy(N,L2_c2(:,3),'k','LineWidth',1,'Marker','*','MarkerSize',6)
grid on
xlabel('M')
ylabel('RMSE')
ylim([1e-10 1e-1])
xlim([4 24])
xticks(4:4:24)
yticks([1e-10 1e-5 1e-1])

nexttile
semilogy(N,nc_c2(:,1),'k','LineWidth',1,'Marker','s','MarkerSize',6)
hold on
semilogy(N,nc_c2(:,2),'k','LineWidth',1,'Marker','^','MarkerSize',6)
semilogy(N,nc_c2(:,3),'k','LineWidth',1,'Marker','*','MarkerSize',6)
grid on
xlabel('M')
ylabel('Condition number')
ylim([1e0 1e22])
xlim([4 24])
xticks(4:4:24)
yticks([1e0 1e5 1e10 1e15 1e20])

allAxes = findall(t, 'Type', 'Axes');
for i = 1:length(allAxes)
    set(allAxes(i), 'FontName', 'Liberation Serif', 'FontSize', 10); % Change 'Arial' to your desired font and 14 to your desired size
    % Update title, xlabel, and ylabel specifically
    % set(get(allAxes(i), 'XLabel'), 'FontName', 'Times', 'FontSize', 12);
    % set(get(allAxes(i), 'YLabel'), 'FontName', 'Times', 'FontSize', 12);
end
leg = legend('Gaussian splines', 'IMQ', 'Wendland');
leg.Layout.Tile = 'north';
legend boxoff  
legend('Orientation','horizontal')
nameOut = 'RBF_comparison';
stmp = strcat('plots/', nameOut, '.pdf');
% exportgraphics(gcf,stmp,'Resolution',400)
exportgraphics(gcf,stmp,'ContentType','vector')





















%% fig INFLUENCE OF DATA POINT LOCATION: uniform vs sinusoidal distribution
L2w = load("L2w.dat");
L2g = load("L2g.dat");
N = 4:2:26;
t = tiledlayout(1,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile
semilogy(N,L2w(:,1),'k','LineWidth',1,'Marker','s','MarkerSize',9)
hold on
semilogy(N,L2w(:,2),'k','LineWidth',1,'Marker','*','MarkerSize',9)
grid on
xlabel('M')
ylabel('RMSE')
ylim([1e-9 1e-2])
xlim([4 24])
xticks(4:4:24)
% set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 12);
% ax = get(gca,'XTickLabel');
% set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 12)
%axis tight


nexttile
semilogy(N,L2g(:,1),'k','LineWidth',1,'Marker','s','MarkerSize',9)
hold on
semilogy(N,L2g(:,2),'k','LineWidth',1,'Marker','*','MarkerSize',9)
%axis tight
ylim([1e-9 1e-2])
xlim([4 24])
xticks(4:4:24)
grid on
xlabel('M')
ylabel('RMSE')
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
leg = legend('Uniform distribution', 'Modified distribution');
leg.Layout.Tile = 'north';
legend boxoff    
legend('Orientation','horizontal')
nameOut = 'dataset';
stmp = strcat('plots/', nameOut, '.pdf');
% exportgraphics(gcf,stmp,'Resolution',400)
exportgraphics(gcf,stmp,'ContentType','vector')