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
% base structure to write xml file
strDomain = readstruct('Domains/domain1block.xml');
interfFile = fullfile('Domains','interfaces_1.xml');
strInterf = readstruct(interfFile);
%% INPUT

% number of refinement
nref = 2;
[h,L2,H1] = deal(zeros(nref,1));

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

  %fprintf('Running mesh refinement %i \n',i);

  % run script to get refined mesh
  fname = "domain_"+elem_type+"_"+num2str(i);

%   command = "python Mesh/scripts/domain.py "  + fname...
%     + " " + num2str(N_i_l) + " " + num2str(N_i_r) + " " + elem_type;
%   system(command);

  meshFile = fullfile('Mesh','meshes',fname+".vtk");
  mesh = Mesh();
  mesh.importMesh(meshFile);

  % set up bc file
  nodes = unique(mesh.surfaces(ismember(mesh.surfaceTag,[2 4]),:));
  c = mesh.coordinates(nodes,:);
  vals = arrayfun(@(i) anal(c(i,:)),1:numel(nodes));
  vals = reshape(vals,[],1);
  writeBCfiles('BCs/bc','NodeBC','Dir','Poisson','manufactured_bc',0,0,nodes,vals);

  %clear mesh

  % write mesh to domain file
  strDomain.Name = "Cube_"+integration_type+"_"+num2str(i);
  strDomain.Geometry = fullfile('Mesh','meshes',fname+".vtk");
  domainFile = fullfile('Domains','domain1block.xml');
  writestruct(strDomain,domainFile);

  % write interface to file
  strInterf.Interface(1).Quadrature.typeAttribute = integration_type;
  strInterf.Interface(1).Quadrature.nGPAttribute = nG;
  %strInterf.Interface(1).Print.nameAttribute = "interf_"+integration_type+"_"+num2str(i);
  if strcmp(integration_type,'RBF')
    strInterf.Interface(1).Quadrature.nIntAttribute = nInt;
  end
  writestruct(strInterf,interfFile);
 
  % processing Poisson problem
  domains = buildModel(domainFile);

  domains.simparams.setVerbosity(0);
  domains.getSolver('Poisson').setAnalSolution(anal,f,gradx,grady,gradz);

  [interfaces,domains] = Mortar.buildInterfaces(interfFile,domains);
  % set up analytical solution

  solver = MultidomainFCSolver(domains,interfaces);
  solver.NonLinearLoop();


  pois = getSolver(domains,'Poisson');
  [L2(i),H1(i)] = pois.computeError_v2();
  h(i) = 1/N_i_r;
  if domains.simparams.verbosity>0
    fprintf('Max absolute error is: %1.6e \n',max(abs(pois.state.data.err)));
  end
end

% compute convergence order
L2ord = log(L2(1:end-1)./L2(2:end))./log(h(1:end-1)./h(2:end));
H1ord = log(H1(1:end-1)./H1(2:end))./log(h(1:end-1)./h(2:end));




