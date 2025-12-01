close all;
% clear;
figures_dir = 'Images';

profile on

scriptFullPath = mfilename('fullpath');
scriptDir = fileparts(scriptFullPath);
cd(scriptDir);
warning('off','MATLAB:nearlySingularMatrix');

% shortcut to define a model using a unique xml file
% useful when dealing with many domains
simparams = SimulationParameters(fullfile("Inputs","domain.xml"));
domain = buildModel(fullfile("Inputs","domain.xml"));

% perform a fully coupled simulation
solver = FCSolver(simparams,domain);
[simState] = solver.NonLinearLoop();

% Finalize the print utility
domain.outstate.finalize()

%% --------------------- Post Processing the Results ----------------------
if true
  image_dir = fullfile(pwd,figures_dir);
  if isfolder(image_dir)
    rmdir(image_dir,"s")
    mkdir(figures_dir)
  else
    mkdir(figures_dir)
  end

  pressure = [domain.outstate.results.pressure];
  displacements = [domain.outstate.results.displacements];
  expTime = [domain.outstate.results.time];

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

  timesInd = [1;2;3];
  time_string = "Year  " + expTime(timesInd);
  set(0,'DefaultAxesColorOrder',[0 0 0],...
    'DefaultAxesLineStyleOrder','-|--|:|-.')

  figure('Position', [100, 100, 700, 700])
  if strcmp()
    pressPlot = pressure(vertNod(indNod),:);
    plot(pressPlot,vertNodZ)
    xlabel('Pressure [kPa]')
    ylabel('z (m)')
    legend(time_string)
  elseif isFVTPFABased(domain.model,'Flow')
    pressPlot = pressure(vertEl(indEl),:);
    plot(pressPlot,vertElZ)
    xlabel('Pressure [kPa]')
    ylabel('z (m)')
    legend(time_string)
  end

  set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
  a = get(gca,'XTickLabel');
  set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10, 'XGrid', 'on', 'YGrid', 'on')

  % export figure with quality
  stmp = fullfile(figures_dir,'DeepAcquifer_pressure.png');
  exportgraphics(gcf,stmp,'Resolution',400)

  dispPlot = displacements(3*vertNod(indNod),:);
  figure('Position', [100, 100, 700, 700])
  plot(1000*dispPlot,vertNodZ);
  xlabel('Vertical displacement (mm)')
  ylabel('z (m)')
  % xlim([0 50])
  % ylim([-60 5])
  legend(time_string)

  set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
  a = get(gca,'XTickLabel');
  set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10, 'XGrid', 'on', 'YGrid', 'on')

  % export figure with quality
  stmp = fullfile(figures_dir,'DeepAcquifer_vertDisplacements.png');
  exportgraphics(gcf,stmp,'Resolution',400)
end