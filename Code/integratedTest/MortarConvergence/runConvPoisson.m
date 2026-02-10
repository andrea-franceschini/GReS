% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);
% Change the current directory to the script's directory
cd(scriptDir);

%% Poisson problem with single domain in 3D. Testing new poisson module

% analytical solution
anal = @(X) cos(pi*X(2)).*cos(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
gradx = @(X) cos(pi*X(2)).*cos(pi*X(3)).*(2 - 2*X(1) + pi*cos(pi*X(1)));
grady = @(X) -pi*sin(pi*X(2)).*cos(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
gradz = @(X) -pi*cos(pi*X(2)).*sin(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
h = @(x) -2-3*pi^2*sin(pi*x)-4*pi^2*x+2*pi^2*x.^2;
f = @(X) cos(pi*X(2)).*cos(pi*X(3)).*h(X(1));


%% model commons
% structure to rewrite the xml file and programatically change the input
fileName = 'poissonMortar.xml';
fileStruct = readstruct(fileName,AttributeSuffix="");



%% INPUT

% number of refinement
nref = 3;
[h,L2,H1] = deal(zeros(nref,1));

% set mortar integration info
nG = 6;
if strcmp(integration_type,'SegmentBased')
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
  meshName = "domain_"+elem_type+"_"+num2str(i);

  % update the mesh in the domain input file
  fileStruct.Domain.Geometry.fileName = fullfile('Input','Mesh','meshes',meshName+".vtk");


  % write interface to file
  fileStruct.Interface.MeshTying.Quadrature.type = integration_type;
  fileStruct.Interface.MeshTying.Quadrature.nGP = nG;
  %strInterf.Interface(1).Print.nameAttribute = "interf_"+integration_type+"_"+num2str(i);
  if strcmp(integration_type,'RBFquadrature')
    fileStruct.Interface(1).MeshTying.Quadrature.nInt = nInt;
  end

  writestruct(fileStruct,fileName,AttributeSuffix="");

  simparams = SimulationParameters(fileName);

  % processing Poisson problem
  [domain,interfaces] = buildModel(fileName);

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




