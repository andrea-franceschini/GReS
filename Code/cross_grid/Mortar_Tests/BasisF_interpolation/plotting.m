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

%% BASIS Function sketch 
% t = tiledlayout(1,2);
% %t.TileSpacing = 'compact';
% t.Padding = 'compact';
% nexttile
figure(1)

x = linspace(-1,1,50);
bf1 = @(x) 0.5*x.*(1+x);
plot(x,bf1(x),'b-','LineWidth',2);
%colormap
xs = linspace(-1,1,5);
hold on
plot([-1 0 1],[0 0 0],'k.-','Linewidth',2,'MarkerSize',30)
scatter(xs,zeros(length(xs),1),"red")
scatter(xs,bf1(xs),"filled","red")
% for i=1:length(xs)
%     plot([xs(i);xs(i)],[0 bf1(xs(i))],'r-','LineWidth',0.5)
% end
xlim([-1.2 1.2])
ylim([-0.3 1.1])
set(gca,'XTick',[])
set(gca,'YTick',[])

%str = {'-1','0','+1'};
%text([-1.15 0 1.05],-0.05*ones(3,1),str,'FontName','Liberation Serif','FontSize', 12)

figure(2)
[X,Y] = meshgrid(-1:.01:1);
bf2 = @(x,y) 0.5*(1-x.^2).*(1+y);
Z = bf2(X,Y);
s = surf(X,Y,Z);
s.EdgeColor = 'none';
s.FaceLighting = "gouraud";
hold on
[xs,ys] = meshgrid(linspace(-1,1,6));
%scatter3([xs(i) xs(i)],[ys(i) ys(i)],[-0.4 bf2(xs(i),ys(i))],"filled","red")
coordLoc = [-1 -1;
   0 -1;
    1 -1;
    1 0;
    1 1;
    0 1;
    -1 1;
    -1 0
    -1 -1];
xp = xs(:);
yp = ys(:);
% for i=1:length(xp)
%     plot3([xp(i) xp(i)],[yp(i) yp(i)],[-0.4 bf2(xp(i),yp(i))],'r-','LineWidth',0.25)
% end
plot3(coordLoc(:,1),coordLoc(:,2),-0.4*ones(length(coordLoc)),'.k-','LineWidth',2,'MarkerSize',30)
scatter3(xs,ys,-0.4*ones(size(xs,1)),"red")
scatter3(xs,ys,bf2(xs,ys),"filled","red")
%str = {'(-1,-1)';'(0,-1)';'(+1,-1)';'(1,0)';'(+1,+1)';'(0,+1)';'(-1,+1)';'(-1,0)'};
cLoc = [-1.2 -1;
   0 -1 ;
    1 -1;
    1.1 0;
    1.1 1.1;
    0 1.4;
    -1.2 1;
    -1.2 0];
%text(cLoc(:,1),cLoc(:,2),-0.45*ones(length(cLoc),1),str,'FontName','Liberation Serif','FontSize', 12)
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'ZTick',[])
xlim([-1.4 1.4])
ylim([-1.1 1.1])
zlim([-0.6 1.2])

allAxes = findall(t, 'Type', 'Axes');
for i = 1:length(allAxes)
    set(allAxes(i), 'FontName', 'Liberation Serif', 'FontSize', 10); % Change 'Arial' to your desired font and 14 to your desired size
    % Update title, xlabel, and ylabel specifically
    % set(get(allAxes(i), 'XLabel'), 'FontName', 'Times', 'FontSize', 12);
    % set(get(allAxes(i), 'YLabel'), 'FontName', 'Times', 'FontSize', 12);
end
stmp = strcat('plots/', 'basisF_sketch', '.png');
% exportgraphics(gcf,stmp,'Resolution',400)
exportgraphics(gcf,stmp,'Resolution',600,'ContentType','image')





















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

%% INFLUENCE OF THE RADIUS
markS = 9;
fix1G = load("graph_rad/fix1g.dat");
% fix05G = load("graph_rad/fix05g.dat");
% fix2G = load("graph_rad/fix2g.dat");
% N = 4:2:24;
% t = tiledlayout(1,2);
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
% nexttile
% semilogy(N,fix05G,'k','LineWidth',1.2,'Marker','s','MarkerSize',markS)
% hold on
% semilogy(N,fix1G,'k','LineWidth',1.2,'Marker','o','MarkerSize',markS)
% semilogy(N,fix2G,'k','LineWidth',1.2,'Marker','^','MarkerSize',markS)
% grid on
% xlabel('M')
% ylabel('RMSE')
% ylim([1e-10 1e-2])
% xlim([4 16])
% xticks(4:4:16)
% legend('r = 0.5 \cdot h_M','r = h_M','r = 2 \cdot h_M')
% set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 12);
% ax = get(gca,'XTickLabel');
% set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 12)
%axis tight

fill2G = load("graph_rad/fill2G.dat");
fill4G = load("graph_rad/fill4G.dat");
fill8G = load("graph_rad/fill8G.dat");
semilogy(N,fill2G,'k--','LineWidth',1.2,'Marker','s','MarkerSize',markS)
hold on
semilogy(N,fill4G,'k--','LineWidth',1.2,'Marker','o','MarkerSize',markS)
semilogy(N,fill8G,'k--','LineWidth',1.2,'Marker','^','MarkerSize',markS)
semilogy(N,fix1G,'k-','LineWidth',1.2,'Marker','*','MarkerSize',markS)
legend('r = 2 \cdot h_{\Xi}','r = 4 \cdot h_{\Xi}','r = 8 \cdot h_{\Xi}','r = h_{M}')
%axis tight
ylim([1e-10 5e-2])
xlim([4 16])
xticks(4:4:16)
grid on
xlabel('M')
ylabel('RMSE')
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 12);
ax = get(gca,'XTickLabel');
set(gca,'XTickLabel',ax,'FontName', 'Liberation Serif','FontSize', 12)

% Set the font name and size for all text in the tiled layout
% allAxes = findall(t, 'Type', 'Axes');
% for i = 1:length(allAxes)
%     set(allAxes(i), 'FontName', 'Liberation Serif', 'FontSize', 12); % Change 'Arial' to your desired font and 14 to your desired size
%     % Update title, xlabel, and ylabel specifically
%     % set(get(allAxes(i), 'XLabel'), 'FontName', 'Times', 'FontSize', 12);
%     % set(get(allAxes(i), 'YLabel'), 'FontName', 'Times', 'FontSize', 12);
% end
% leg = legend('Uniform distribution', 'Modified distribution');
% leg.Layout.Tile = 'north';
% legend boxoff    
% legend('Orientation','horizontal')
nameOut = 'radius';
stmp = strcat('plots/', nameOut, '.pdf');
% exportgraphics(gcf,stmp,'Resolution',400)
exportgraphics(gcf,stmp,'ContentType','vector')

