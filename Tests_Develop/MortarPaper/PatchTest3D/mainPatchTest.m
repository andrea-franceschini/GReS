% simple patch test to validate hexa 27 finite elements

close all;
clear;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
cd(scriptDir)

integration_type = ["ElementBased","SegmentBased","RBF"];
%nG = [4,8,16];
nInt = 2;

outStr = [];
% write interface to file

for i_t = integration_type

  if strcmp(i_t,'SegmentBased')
    nG = 7;
  else
    nG = [4 16];
  end
  for ng = nG
    interfFile = 'interface.xml';
    strInterf = readstruct(interfFile);
    strInterf.Interface(1).Quadrature.typeAttribute = i_t;
    strInterf.Interface(1).Quadrature.nGPAttribute = ng;
    strInterf.Interface(1).Print.nameAttribute = "interf_"+i_t;
    if strcmp(i_t,'RBF')
      strInterf.Interface(1).Quadrature.nIntAttribute = nInt;
    end

    writestruct(strInterf,interfFile);

    % Set physical models
    model = ModelType("Poromechanics_FEM");

    % Set parameters of the simulation
    fileName = "simParam.dat";
    simParam = SimulationParameters(fileName,model);

    top = Mesh();
    top.importMesh(fullfile('Mesh','top.vtk'));

    bot = Mesh();
    bot.importMesh(fullfile('Mesh','bottom.vtk'));

    writeBCfiles('BCs/bc_bot','SurfBC','Dir',{'Poromechanics','x','y','z'},'bottom_fixed',0,0,bot,1);
    writeBCfiles('BCs/bc_top','SurfBC','Neu',{'Poromechanics','z'},'top_load',0,-0.5,top,2); % left block lateral fix

    % processing Poisson problem
    domains = buildModelStruct_new('domain2block.xml',simParam);

    [interfaces,domains] = Mortar.buildInterfaceStruct('interface.xml',domains);
    % set up analytical solution

    solver = MultidomainFCSolver(simParam,domains,interfaces);
    solver.NonLinearLoop();
    solver.finalizeOutput();


    % collect multipliers on the diagonal of the slave domain
    mult = solver.interfaces{1}.multipliers(1).curr(3*[1 13 16 3]);

    outStr = [outStr;...
      struct('integration',i_t,'nG',ng,'mult',mult)];
  end
end


%%
figure(1)
clf
diagcoord = sqrt(2)*[0 1/3 2/3 1];
hold on

legendEntries = {};  % per raccogliere le voci della legenda

for i = 1:numel(outStr)
  switch outStr(i).integration
    case 'RBF'
      lc = 'b';
    case 'ElementBased'
      lc = 'g';
    case 'SegmentBased'
      lc = 'r';
    otherwise
      lc = 'k';  % colore di fallback
  end
  mult = outStr(i).mult;
  h = plot(diagcoord, mult, '-s', 'Color', lc, 'LineWidth', 1.5, ...
           'DisplayName', outStr(i).integration);
end

% Etichette assi in LaTeX
xlabel('coordinate diagonal', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 14)

% Legenda
%legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'best')

% Font LaTeX per tutto il grafico
set(gca, 'FontName', 'latex', 'TickLabelInterpreter', 'latex')
grid on
box on



