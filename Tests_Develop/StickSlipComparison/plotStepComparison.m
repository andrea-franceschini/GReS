%% Plot tractions along crack for requested time step
tStep = 3;

% --- Files and labels
historyFiles = {...
    fullfile("OUT/historyCrackP0.mat"),...
    fullfile("OUT/historyCrackBubble.mat"),...
    fullfile("OUT/historyCrackDual.mat"),...
    fullfile("OUT/historyCrackStandard.mat"),...
    fullfile("OUT/historyCrackFine.mat")...
    };

legendLabels = {...
    'P0 stabilized', ...
    'Bubble stabilized', ...
    'Dual multipliers', ...
    'Standard multipliers', ...
    'Dual conforming (reference)'...
    };

% --- Colors (RGB)
colors = [ ...
    1.0, 0.000, 0.0;   % magenta (P0 stabilized)  <-- MAIN
    0.850, 0.550, 0.000;   % muted orange (bubble)
    0.300, 0.650, 0.300;   % muted green (dual)
    0.0, 0.0, 1.0;   % muted red (standard)
    0.000, 0.000, 0.000;   % black (fine conforming)
    ];

% --- Line styles
lineStyles = {'-', '--', ':', '-.', '-'};

% --- Markers
markers = {'o','s','^','v','none'};

% --- Line widths
lineWidths  = [2.2, 1.0, 1.0, 1.0, 2.2];   % emphasize P0 + fine
markerSizes = [6, 3.5, 3.5, 3.5, 0.01];       % secondary = less invasive

% --- Initialize figures
figure(1); clf; hold on;
xlabel('$\sigma_n$', 'Interpreter','latex','FontSize',14);
ylabel('$z$ (m)', 'Interpreter','latex','FontSize',14);

figure(2); clf; hold on;
xlabel('Norm tangential traction', 'Interpreter','latex','FontSize',14);
ylabel('$z$ (m)', 'Interpreter','latex','FontSize',14);

% --- Loop over files
for k = 1:length(historyFiles)
    historyFile = historyFiles{k};
    vars = load(historyFile);
    mult = [vars.output.multipliers];
    
    mult = mult(:,tStep);
    sn = mult(1:3:end);
    norm_tT = sqrt(mult(2:3:end).^2 + mult(3:3:end).^2);
    
    interf = vars.interfaces{1};
    mshSlave = getMesh(interf, MortarSide.slave);
    surfC = mshSlave.surfaceCentroid;
    
    NX = round(10/(2*surfC(1,2)));
    NY = round(15/(2*surfC(1,3)));
    
    if interf.multiplierLocation == entityField.surface
        vertAxis1 = 0.5*NX:NX:(NX*NY - 0.5*NX);
        vertAxis2 = vertAxis1+1;
        
        snPlot = 0.5*(sn(vertAxis1)+sn(vertAxis2));
        tTPlot = 0.5*(norm_tT(vertAxis1)+norm_tT(vertAxis2));
        zPlot = surfC(vertAxis2,3);
    else
        [~,i] = sort(mshSlave.coordinates(:,3));
        coord = mshSlave.coordinates(i,:);
        nodes = coord(:,2) == 5;
        sn = sn(i);
        tT = norm_tT(i);
        snPlot = sn(nodes);
        tTPlot = tT(nodes);
        zPlot = linspace(0,15,NY+1);
    end
    
    % --- Plot normal traction
    figure(1);
    plot(snPlot, zPlot, ...
        'Color', colors(k,:), ...
        'LineStyle', lineStyles{k}, ...
        'Marker', markers{k}, ...
        'LineWidth', lineWidths(k), ...
        'MarkerSize', markerSizes(k), ...
        'MarkerFaceColor', colors(k,:));
    
    % --- Plot tangential traction
    figure(2);
    plot(tTPlot, zPlot, ...
        'Color', colors(k,:), ...
        'LineStyle', lineStyles{k}, ...
        'Marker', markers{k}, ...
        'LineWidth', lineWidths(k), ...
        'MarkerSize', markerSizes(k), ...
        'MarkerFaceColor', colors(k,:));
end

% --- Set axes ticks, fonts, grid
figs = [1,2];
for f = figs
    figure(f);
    set(gca, 'TickLabelInterpreter','latex', 'FontSize',12, 'LineWidth',1.1);
    grid on; box on;
end

% --- Add legends
figure(1); legend(legendLabels,'Interpreter','latex','Location','best');
figure(2); legend(legendLabels,'Interpreter','latex','Location','best');

exportgraphics(figure(1), 'OUT/tractions_normal.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');

exportgraphics(figure(2), 'OUT/tractions_tangential.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');