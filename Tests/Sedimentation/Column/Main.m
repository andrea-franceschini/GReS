% close all;
% clear;
input_dir = 'Inputs/';
file_Mat = fullfile(input_dir,'materialsB.xml');
% profile on;

solv = {'A','B','C','D'};

%% ------------------------------ Set up the Domain -----------------------
% Set the simulation parameters for the non-linear solver.
simParam = SimulationParameters('Start',0.,'End',400.,...
      'DtInit',1,'DtMin',1e-2,'DtMax',1);

% Create an object of the Materials class and read the materials file
mat = Materials(file_Mat);

cells2Plot = [4];
sim = struct;
for mod = 1:numel(solv)
  file_Solver = fullfile(input_dir,"solver"+solv(mod)+".xml");
  sol = fullfile(input_dir,"solver"+solv(mod)+".xml");

  printUtils = OutState('outputFile',strcat('Outputs/Results', solv{mod}),'printTimes',3000);
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
% profile off;
% profile viewer;


%% -------------------------- Print some comparison -----------------------
tstr = "Case " + string(solv(:));

if false
  figure('Position',[100 100 700 700])
  hold on
  for mod = 1:numel(solv)
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
  for mod = 1:numel(solv)
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
  for mod = 1:numel(solv)
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