close all;
% clear;
input_dir = 'Inputs/';
file_SimP = fullfile(input_dir,'simparam.xml');
file_Mat = fullfile(input_dir,'materials.xml');
file_Output = fullfile(input_dir,'output.xml');
% profile on;

solver_files = {'solverA.xml','solverB.xml','solverC.xml','solverD.xml'};
% solver_files = {'solverD.xml'};

%% ------------------------------ Set up the Domain -----------------------
% Set the simulation parameters for the non-linear solver.
simParam = SimulationParameters(file_SimP);

% Create an object of the Materials class and read the materials file
mat = Materials(file_Mat);

cells2Plot = [4];
sim = struct;
for mod = 1:numel(solver_files)  
  file_Solver = fullfile(input_dir, solver_files{mod});

  printUtils = OutState(file_Output);
  helpFunc = @(f) cell2mat(arrayfun(@(s) s.(f)(cells2Plot), ...
                      printUtils.results, 'UniformOutput', false))';
  
  % Create object handling construction of Jacobian and rhs of the model
  domain = Discretizer('Materials',mat);
  domain.addPhysicsSolvers(file_Solver);

  % The modular structure of the discretizer class allow the user to easily
  % customize the solution scheme.
  solver = EvolvingGrid('simulationparameters',simParam,...
    'domains',domain,...
    'output',printUtils, ...
    'growprint',0, ... 
    'intervalprint',[0,400]);

  solver.simulationLoop();

  sim(mod).pressure = helpFunc('pressure');
  sim(mod).stress   = helpFunc('stress');
  sim(mod).strain   = helpFunc('strain');
  sim(mod).porosity = helpFunc('porosity');
  sim(mod).t = [printUtils.results.time];
end

%%
tstr = compose("Case %d", 1:numel(solver_files))';

if true

figure('Position',[100 100 700 700])
hold on
for mod = 1:numel(solver_files)  
  plot(sim(mod).t, sim(mod).pressure, '+-', 'LineWidth', 2);
  tmin = floor(min(sim(mod).t));
  tmax = ceil(max(sim(mod).t));
end

xlabel('Time (year)')
ylabel('Pressure (Pa)')
legend(tstr, 'Location','northwest')

ax = gca;
set(ax, 'FontName','Liberation Serif', 'FontSize',16, ...
  'XLim',[tmin tmax], 'XTick',0:5:tmax, 'XMinorTick','on', ...
  'XGrid','on', 'YGrid','on')
ax.XAxis.MinorTickValues = 0:1:tmax;
grid minor

end

if false

figure('Position',[100 100 700 700])
hold on
for mod = 1:numel(solver_files)  
  plot(sim(mod).t, -sim(mod).stress, '+-', 'LineWidth', 2);
  tmin = floor(min(sim(mod).t));
  tmax = ceil(max(sim(mod).t));
end

ylabel('Stress (Pa)')

xlabel('Time (year)')
legend(tstr, 'Location','northwest')

ax = gca;
set(ax, 'FontName','Liberation Serif', 'FontSize',16, ...
  'XLim',[tmin tmax], 'XTick',0:5:tmax, 'XMinorTick','on', ...
  'XGrid','on', 'YGrid','on')
ax.XAxis.MinorTickValues = 0:1:tmax;
grid minor

end



if false

figure('Position',[100 100 700 700])
hold on
for mod = 1:numel(solver_files)  
  plot(sim(mod).t, 100*sim(mod).porosity, '+-', 'LineWidth', 2);
  tmin = floor(min(sim(mod).t));
  tmax = ceil(max(sim(mod).t));
end

ylabel('Porosity (%)')

xlabel('Time (year)')
legend(tstr, 'Location','northwest')

ax = gca;
set(ax, 'FontName','Liberation Serif', 'FontSize',16, ...
  'XLim',[tmin tmax], 'XTick',0:5:tmax, 'XMinorTick','on', ...
  'XGrid','on', 'YGrid','on')
ax.XAxis.MinorTickValues = 0:1:tmax;
grid minor

end













% profile off;
% profile viewer;
%%
postproc=false;
if postproc
  cells2Plot = [1;4];
  outFolder = "Outputs/Figs";
  image_dir = fullfile(pwd,outFolder);
  if isfolder(image_dir)
    rmdir(image_dir,"s")
    mkdir(outFolder)
  else
    mkdir(outFolder)
  end

  % Saving a temporary variabel.
  helpFunc = @(f) cell2mat(arrayfun(@(s) s.(f)(cells2Plot), ...
                      printUtils.results, 'UniformOutput', false))';
  pressure = helpFunc('pressure');
  stress   = helpFunc('stress');
  strain   = helpFunc('strain');
  porosity = helpFunc('porosity');
  t = [printUtils.results.time];
  % save('dh40.mat', 'pressure', 'stress', 'strain', 'porosity', 't');

  % Ajusting the legend.
  tstr = compose("Point %d", 1:length(cells2Plot))';

  if true
    %Plotting pressure head
    figure('Position',[100 100 700 700])
    plot(t, pressure, '+-', 'LineWidth', 2);
    hold on
    xlabel('Time (year)')
    ylabel('Pressure (Pa)')
    legend(tstr, 'Location','northwest')

    tmin = floor(min(t));
    tmax = ceil(max(t));
    ax = gca;
    set(ax, 'FontName','Liberation Serif', 'FontSize',16, ...
      'XLim',[tmin tmax], 'XTick',0:5:tmax, 'XMinorTick','on', ...
      'XGrid','on', 'YGrid','on')
    ax.XAxis.MinorTickValues = 0:1:tmax;
    grid minor

    % export figure with quality
    stmp = fullfile(image_dir,'pressure.png');
    exportgraphics(gcf,stmp,'Resolution',400)
  end

  if false
    figure('Position', [100, 100, 700, 700])
    plot(t,-stress,'-', 'LineWidth', 2, 'MarkerSize', 14);
    hold on
    ylabel('Stress (Pa)')
    xlabel('Time (year)')
    legend(tstr, 'Location', 'northwest')
    set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
    % export figure with quality
    stmp = fullfile(image_dir,'stress.png');
    exportgraphics(gcf,stmp,'Resolution',400)
  end

  if false
    figure('Position', [100, 100, 700, 700])
    plot(t,-100*strain,'-', 'LineWidth', 2, 'MarkerSize', 14);
    hold on
    ylabel('Strain (%)')
    xlabel('Time (year)')
    legend(tstr, 'Location', 'northwest')
    set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
    % export figure with quality
    stmp = fullfile(image_dir,'strain.png');
    exportgraphics(gcf,stmp,'Resolution',400)
  end

  if false
    figure('Position', [100, 100, 700, 700])
    plot(t,porosity,'-', 'LineWidth', 2, 'MarkerSize', 14);
    hold on
    ylabel('Porosity (-)')
    xlabel('Time (year)')
    legend(tstr, 'Location', 'northeast')
    set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
    % export figure with quality
    stmp = fullfile(image_dir,'porosity.png');
    exportgraphics(gcf,stmp,'Resolution',400)
  end
end




if false
  load('dh.mat');
  tstr = compose("Cell Height %.1f", [0.05 0.1 0.2 0.4])';
  listVar=[dh05,dh10,dh20,dh40];
  figure('Position', [100, 100, 700, 700])
  hold on
  for sim=1:length(listVar)
    plot(listVar(sim).t,listVar(sim).pressure(:,2),'-', 'LineWidth', 2, 'MarkerSize', 14);
  end
  ylabel('Pressure (Pa)')
  xlabel('Time (year)')
  legend(tstr, 'Location', 'northwest')
  set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
end