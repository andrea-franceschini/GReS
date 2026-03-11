% study parameters

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);
% Change the current directory to the script's directory
cd(scriptDir);

fileName = 'poissonMortar.xml';
params = readstruct(fileName,AttributeSuffix="");

for elem_type = ["hexa","hexa27"]

  for integration_type = ["SegmentBasedQuadrature",...
                          "RBFquadrature",...
                          "ElementBasedQuadrature",...
                          ]

    [L2,H1] = run(params,elem_type,integration_type);
    validate(elem_type,L2,H1);

    clearvars

  end
end


function [L2ord,H1ord] = run(params,elem,quadrature)


%% Poisson problem with single domain in 3D. Testing new poisson module

% analytical solution
anal = @(X) cos(pi*X(2)).*cos(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
gradx = @(X) cos(pi*X(2)).*cos(pi*X(3)).*(2 - 2*X(1) + pi*cos(pi*X(1)));
grady = @(X) -pi*sin(pi*X(2)).*cos(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
gradz = @(X) -pi*cos(pi*X(2)).*sin(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
h = @(x) -2-3*pi^2*sin(pi*x)-4*pi^2*x+2*pi^2*x.^2;
f = @(X) cos(pi*X(2)).*cos(pi*X(3)).*h(X(1));



%% INPUT

% number of refinement
nref = 2;
[h,L2,H1] = deal(zeros(nref,1));

% set mortar integration info
nG = 6;
if strcmp(quadrature,'SegmentBased')
  nG = 7;
end

nInt = 5;

N_l = [4 8 16];

N_r = [6 12 24];

%% convergence loop
for i = 1:nref
  N_i_l = N_l(i);
  N_i_r = N_r(i);

  % run script to get refined mesh
  meshName = "domain_"+elem+"_"+num2str(i);
  meshFileName = fullfile('Input','Mesh','meshes',meshName+".vtk");

  mesh = Mesh();
  mesh.importMesh(meshFileName);
  elems = Elements(mesh,3);
  grid = struct('topology',mesh,'cells',elems);
  bc = Boundaries(grid);
  bc.addBC('name',"ManufacturedSolution")


  domain 


  
  % write interface to file
  params.Interface.MeshTying.Quadrature.type = quadrature;
  params.Interface.MeshTying.Quadrature.nGP = nG;
  %strInterf.Interface(1).Print.nameAttribute = "interf_"+integration_type+"_"+num2str(i);
  if strcmp(quadrature,'RBFquadrature')
    params.Interface(1).MeshTying.Quadrature.nInt = nInt;
  end

  simparams = SimulationParameters(params.SimulationParameters);

  % processing Poisson problem

  [domain,interfaces] = buildModel(params);

  domain.getPhysicsSolver("Poisson").setAnalSolution(anal,f,gradx,grady,gradz);

  solver = NonLinearImplicit('simulationparameters',simparams,...
                             'domains',domain,...
                             'interface',interfaces);
  solver.simulationLoop();

  pois = getPhysicsSolver(domain,'Poisson');
  [L2(i),H1(i)] = pois.computeError_v2();
  h(i) = 1/N_i_r;

  gresLog().log(0,'Max absolute error is: %1.6e \n',max(abs(domain.state.data.err)));
end

% compute convergence order
L2ord = log(L2(1:end-1)./L2(2:end))./log(h(1:end-1)./h(2:end));
H1ord = log(H1(1:end-1)./H1(2:end))./log(h(1:end-1)./h(2:end));



end

function validate(elem,L2,H1)
    msg = "Error for %s element with %s mortar quadrature scheme";

    switch elem
      % ensure correct convergence rate 
      case "hexa"
        assert(all([L2>1.8;L2<2.5]),...
          msg,elem_type,integration_type)
        assert(all([H1>0.8;H1<1.5]),...
          msg,elem_type,integration_type)
      case "hexa27"
        assert(all([L2>2.8;L2<3.5]),...
          msg,elem_type,integration_type)
        assert(all([H1>1.8;H1<2.5]),...
          msg,elem_type,integration_type)
    end
end
