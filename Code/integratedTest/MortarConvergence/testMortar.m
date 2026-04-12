% study parameters

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);
% Change the current directory to the script's directory
cd(scriptDir);


for elem_type = ["hexa","hexa27"]

  for integration_type = ["SegmentBasedQuadrature",...
                          "ElementBasedQuadrature",...
                          "RBFquadrature"
                          ]

    [L2,H1] = run(elem_type,integration_type);
    validate(elem_type,integration_type,L2,H1);

  end
end


function [L2ord,H1ord] = run(elem,quadrature)

fileName = 'poissonMortar.xml';
params = readstruct(fileName,AttributeSuffix="");


%% Poisson problem with single domain in 3D. Testing new poisson module

% analytical solution
anal = @(x,y,z) cos(pi*y).*cos(pi*z).*(2*x-x.^2 + sin(pi*x));
gradx = @(x,y,z) cos(pi*y).*cos(pi*z).*(2 - 2*x + pi*cos(pi*x));
grady = @(x,y,z) -pi*sin(pi*y).*cos(pi*z).*(2*x-x.^2 + sin(pi*x));
gradz = @(x,y,z) -pi*cos(pi*y).*sin(pi*z).*(2*x-x.^2 + sin(pi*x));
h = @(x,y,z) -2-3*pi^2*sin(pi*x)-4*pi^2*x+2*pi^2*x.^2;
f = @(x,y,z) cos(pi*y).*cos(pi*z).*h(x);



%% INPUT

% number of refinement
nref = 2;
[h,L2,H1] = deal(zeros(nref,1));

% set mortar integration info
if strcmp(quadrature,'SegmentBasedQuadrature')
  nGP = 7;
else
  nGP = 6;
end

interfStr.masterDomain = 1;
interfStr.slaveDomain = 1;
interfStr.masterSurface = 1;
interfStr.slaveSurface = 3;
interfStr.multiplierType="dual";

quadStr.type = quadrature;
quadStr.nGP = nGP;
quadStr.nInt = 5;


%Nl = [4 8 16];
Nr = [6 12 24];

simparams = SimulationParameters(params.SimulationParameters);


%% convergence loop
for i = 1:nref
  %N_i_l = N_l(i);

  % run script to get refined mesh
  meshName = "domain_"+elem+"_"+num2str(i);
  meshFileName = fullfile('Input','Mesh','meshes',meshName+".vtk");

  grid = Grid();
  grid.importMesh(meshFileName);

 surf = grid.surfaces;


  bcEnts = unique(surf.connectivity(ismember(surf.tag,[2;4]),:));
  c = grid.coordinates(bcEnts,:);
  [X,Y,Z] = deal(c(:,1),c(:,2),c(:,3)); 
  bcVals = anal(X,Y,Z);
  bc = Boundaries(grid);
  bc.addBC('name',"analSol",...
    'type',"dirichlet",...
    'field',"node",...
    'entityListType',"bcList",...
    'entityList',bcEnts,...
    'variable',"u")

  bc.addBCEvent("analSol",'time',0.0,'value',bcVals)

  domain = Discretizer('grid',grid,'boundaries',bc);
  domain.addPhysicsSolver('Poisson');

  interfStr.Quadrature = quadStr;

  interface = InterfaceSolver.add('MeshTying',domain,interfStr);

  domain.getPhysicsSolver("Poisson").setAnalSolution(anal,f,gradx,grady,gradz);

  solver = NonLinearImplicit('simulationparameters',simparams,...
    'domains',domain,...
    'interface',interface);

  solver.simulationLoop();

  pois = getPhysicsSolver(domain,'Poisson');
  [L2(i),H1(i)] = pois.computeError_v2();
  h(i) = 1/Nr(i);

  gresLog().log(0,'Max absolute error is: %1.6e \n',max(abs(domain.state.data.err)));
end

% compute convergence order
L2ord = log(L2(1:end-1)./L2(2:end))./log(h(1:end-1)./h(2:end));
H1ord = log(H1(1:end-1)./H1(2:end))./log(h(1:end-1)./h(2:end));



end

function validate(elem,quadrature,L2,H1)
    msg = "Error for %s element with %s mortar quadrature scheme";

    switch elem
      % ensure correct convergence rate 
      case "hexa"
        assert(all([L2>1.8;L2<2.5]),...
          msg,elem,quadrature)
        assert(all([H1>0.8;H1<1.5]),...
          msg,elem,quadrature)
      case "hexa27"
        assert(all([L2>2.8;L2<3.5]),...
          msg,elem,quadrature)
        assert(all([H1>1.8;H1<2.5]),...
          msg,elem,quadrature)
    end
end
