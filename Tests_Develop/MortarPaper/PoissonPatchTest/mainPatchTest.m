% simple patch test to validate hexa 27 finite elements

close all;
%clear;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
cd(scriptDir)

integration_type = ["SegmentBased"];
%nG = [4,8,16];
nInt = 4;

% analytical solution
anal = @(X) X(3);
gradx = @(X) 0;
grady = @(X) 0;
gradz = @(X) 1;
f = @(X) 0;

outStr = [];
% write interface to file

for i_t = integration_type

  if strcmp(i_t,'SegmentBased')
    nG = 7;
  else
    nG = [4];
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
    model = ModelType("Poisson_FEM");

    % Set parameters of the simulation
    fileName = "simParam.dat";
    simParam = SimulationParameters(fileName,model);

    top = Mesh();
    top.importMesh(fullfile('Mesh','top.msh'));

    bot = Mesh();
    bot.importMesh(fullfile('Mesh','bottom.msh'));


    % set up bc file top
    nodes = unique(top.surfaces(ismember(top.surfaceTag,1),:));
    c = top.coordinates(nodes,:);
    vals = arrayfun(@(i) anal(c(i,:)),1:numel(nodes));
    vals = reshape(vals,[],1);
    writeBCfiles('BCs/bc_top','NodeBC','Dir','Poisson','manufactured_solution_top',0,0,nodes,vals);

    % set up bc file bottom
    nodes = unique(bot.surfaces(ismember(bot.surfaceTag,1),:));
    c = bot.coordinates(nodes,:);
    vals = arrayfun(@(i) anal(c(i,:)),1:numel(nodes));
    vals = reshape(vals,[],1);
    writeBCfiles('BCs/bc_bot','NodeBC','Dir','Poisson','manufactured_solution_bot',0,0,nodes,vals);
 
    % processing Poisson problem
    domains = buildModelStruct_new('domain2block.xml',simParam);
    domains(1).Discretizer.getSolver('Poisson').setAnalSolution(anal,f,gradx,grady,gradz);
    domains(2).Discretizer.getSolver('Poisson').setAnalSolution(anal,f,gradx,grady,gradz);

    [interfaces,domains] = Mortar.buildInterfaceStruct('interface.xml',domains);
    % set up analytical solution

    solver = MultidomainFCSolver(simParam,domains,interfaces);
    solver.NonLinearLoop();
    solver.finalizeOutput();


    % collect multipliers on the diagonal of the slave domain
    mult = solver.interfaces{1}.multipliers(1).curr([1 13 16 3]);

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
      mark = 's';
    case 'ElementBased'
      lc = 'r';
      mark = '^';
    case 'SegmentBased'
      lc = 'k';
      mark = 'o';
    otherwise
      lc = 'k';  % colore di fallback
  end
  mult = outStr(i).mult;
  h = plot(diagcoord, mult, 'Color', lc, 'LineWidth', 1.2,'Marker',mark, ...
           'Markersize',9,'LineStyle','-','DisplayName', outStr(i).integration);
end

% Etichette assi in LaTeX
xlabel('coordinate diagonal', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 14)

% Legenda
legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'best')

% Font LaTeX per tutto il grafico
set(gca, 'FontName', 'latex', 'TickLabelInterpreter', 'latex')
grid on
box on
exportgraphics(gcf, 'patchTest.pdf', 'ContentType', 'vector');



