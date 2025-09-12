close all;
% clear;

profile on

scriptFullPath = mfilename('fullpath');
scriptDir = fileparts(scriptFullPath);
cd(scriptDir);
warning('off','MATLAB:nearlySingularMatrix');

% shortcut to define a model using a unique xml file
% useful when dealing with many domains
domain = buildModel('domain.xml');

% perform a fully coupled simulation
solver = FCSolver(domain);
[simState] = solver.NonLinearLoop();

% Finalize the print utility
domain.outstate.finalize()

%% POST PROCESSING

image_dir = fullfile(pwd,'Images');
if isfolder(image_dir)
    rmdir(image_dir,"s")
    mkdir Images
else
    mkdir Images
end

expPress = [domain.outstate.results.expPress];
expDispl = [domain.outstate.results.expDispl];
expTime = [domain.outstate.results.expTime];

topol = domain.grid.topology;

%find nodes in vertical symmetry axis
tmp1=topol.coordinates(:,1)<500.1;
tmp2 = topol.coordinates(:,1)>499.9;
tmp3 = topol.coordinates(:,2)<500.1;
tmp4 = topol.coordinates(:,2)>499.9;
tmpNod = tmp1+tmp2+tmp3+tmp4;
vertNod = find(tmpNod == 4);
[vertNodZ,indNod] = sort(topol.coordinates(vertNod,3));

%find elemes in vertical symmetry axis
tmp1 = topol.cellCentroid(:,1)<450.1;
tmp2 = topol.cellCentroid(:,1)>449.9;
tmp3 = topol.cellCentroid(:,2)<550.1;
tmp4 = topol.cellCentroid(:,2)>449.9;
tmpEl = tmp1+tmp2+tmp3+tmp4;
vertEl = find(tmpEl == 4);
[vertElZ,indEl] = sort(topol.cellCentroid(vertEl,3));

timesInd = [2;3;4];
time_string = "Year  " + expTime(timesInd);
set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|--|:|-.')

figure('Position', [100, 100, 700, 700])
if isFEMBased(domain.model,'Flow')
    pressPlot = expPress(vertNod(indNod),:);
    plot(pressPlot,vertNodZ)
    xlabel('Pressure [kPa]')
    ylabel('z (m)')
    legend(time_string)
elseif isFVTPFABased(domain.model,'Flow')
    pressPlot = expPress(vertEl(indEl),:);
    plot(pressPlot,vertElZ)
    xlabel('Pressure [kPa]')
    ylabel('z (m)')
    legend(time_string)
end
grid on
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
% export figure with quality
stmp = fullfile(image_dir,'DeepAcquifer_pressure.png');
exportgraphics(gcf,stmp,'Resolution',400)



dispPlot = expDispl(3*vertNod(indNod),:);

figure('Position', [100, 100, 700, 700])
% figure(2)
plot(1000*dispPlot,vertNodZ);
xlabel('Vertical displacement (mm)')
ylabel('z (m)')
% xlim([0 50])
% ylim([-60 5])
legend(time_string)
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
% export figure with quality
stmp = fullfile(image_dir,'DeepAcquifer_vertDisplacements.png');
exportgraphics(gcf,stmp,'Resolution',400)
