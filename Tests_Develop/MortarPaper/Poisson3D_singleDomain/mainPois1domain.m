clear
close all

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
% Set physical models 
model = ModelType("Poisson_FEM");
% Set parameters of the simulation
fileName = "simParam.dat";
simParam = SimulationParameters(fileName,model);
% Create an object of the Materials class and read the materials file
mat = [];
% Create the Mesh object
topology = Mesh();

%% INPUT
% base domain size
N = 3;
% number of refinement
nref = 2;
[h,L2,H1] = deal(zeros(nref,1));

elemType = 'hexa27';

switch elemType
  case 'tetra'
    nG = 1;
  case 'hexa'
    nG = 2;
  case 'hexa27'
    nG = 3;
end


%% convergence loop 
for i = 1:nref
  N_i = N*2^(i-1);
  fname = strcat('domain_',num2str(i));
    command = "python Mesh/SingleDomain.py "  + fname...
    + " " + num2str(N_i) + " " + elemType;
  system(command);
  fprintf('Running mesh refinement %i \n',i);
  runPoisson;

  pois = getSolver(linSyst,'Poisson');
  [L2(i),H1(i)] = pois.computeError_v2();
  h(i) = 1/N_i;
  fprintf('Max absolute error is: %1.6e \n',max(abs(pois.state.data.err)));
end

% compute convergence order
L2ord = log(L2(1:end-1)./L2(2:end))./log(h(1:end-1)./h(2:end));
H1ord = log(H1(1:end-1)./H1(2:end))./log(h(1:end-1)./h(2:end));

%% plotting convergence profiles
figure(1)
loglog(h,L2,'-ro')
hold on
loglog(h,H1,'-b^')
xlabel('h')
ylabel('error_norm')
legend('L2','H1')




