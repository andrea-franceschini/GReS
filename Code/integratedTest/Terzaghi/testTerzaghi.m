clear


% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);
% Change the current directory to the script's directory
cd(scriptDir);

for elem = ["hexa","tetra"]
  for flow = ["FV","FEM"]
    if elem == "tetra" && flow == "FV"
      continue
    end
    for scheme = ["FC","FS"]
      run(elem,flow,scheme);
    end
  end
end


function run(elemShape,flowSolver,scheme)

switch elemShape
  case 'tetra'
    ng = 1;
    topology = Mesh();
    topology.importGMSHmesh('Input/Mesh/Column_tetra.msh');
  case 'hexa'
    ng = 2;
    topology = structuredMesh(1,1,50,[0 1],[0 1],[0 10]);
end

simParam = SimulationParameters("Input/simParam.xml");
mat = Materials('Input/materials.xml');
elems = Elements(topology,ng);

if strcmp(flowSolver,"FV")
  solverIn = struct("SinglePhaseFlowFVTPFA",[]);
    faces = Faces(topology);
  grid = struct('topology',topology,'cells',elems,'faces',faces);
else
  solverIn = struct("SinglePhaseFlowFEM",[]);
  grid = struct('topology',topology,'cells',elems);
end

printUtils = OutState('Input/output.xml');

% Create an object of the "Boundaries" class 
bound = Boundaries(grid,'Input/boundaryConditions.xml');

% create the Discretizer (key-value pair input)
domain = Discretizer('grid',grid,...
                     'materials',mat,...
                     'boundaries',bound);

switch scheme
  case 'FC'
    domain.addPhysicsSolver('BiotFullyCoupled',solverIn);
  case 'FS'
    domain.addPhysicsSolver('BiotFixedStressSplit',solverIn);
end

% manually apply initial conditions
state = domain.getState();
applyTerzaghiIC(state,mat,topology,-10);

switch scheme
  case 'FC'
    solver = NonLinearImplicit('simulationparameters',simParam,...
      'domains',domain,...
      'output',printUtils);
  case 'FS'
    solver = FixedStressSplit('simulationparameters',simParam,...
                              'domains',domain,...
                              'output',printUtils,...
                              'maxiterations',20,...
                              'reltolerance',1e-7);

end

solver.simulationLoop();

validate(solver)

end



function validate(solver)

[pAnalytical,uAnalytical] = computeAnalyticalSolution(solver,-10);

pNumerical = [solver.output.results.pressure];
uNumerical = [solver.output.results.displacements];


for i = 1:numel(solver.output.timeList)
  % compare pressure and displacement solution at each time step

  pAn = abs(pAnalytical(:,i));
  pNum = pNumerical(:,i);

  uAn = uAnalytical(:,i);
  uNum = uNumerical(:,i);
  uNum = uNum(3:3:end);

  relErrP = (pAn-pNum)./pNum;
  relErrP(isinf(relErrP)) = 0;
  relErrP(isnan(relErrP)) = 0;
  assert(norm(relErrP)<1e0,'Pressure error for out time %i \n',i)

  relErrU = (uAn-uNum)./(uNum);
  relErrU(isinf(relErrU)) = 0;
  assert(norm(relErrU)<1e0,'Displacement error for out time %i\n',i)
end

if strcmp(class(solver),"FixedStressSplit")
  nIter = mean(solver.getFixedStressIer);
  assert(nIter==2,"Unexpected number of fixed stress split iterations")
end

end




function [p,u] = computeAnalyticalSolution(solv,F)
%fprintf('Computing Terzaghi Analytical solution... \n');
% number of terms for analyitcal soluton series
nm = 1000;

mat = solv.domains.materials;
grid = solv.domains.grid;
mesh = grid.topology;


% Get model geometry
L_min = min(mesh.coordinates(:,3));
L_max = max(mesh.coordinates(:,3));
L = abs(L_max-L_min);

% Get Material parameters from materials class
k = mat.getMaterial(1).PorousRock.getPermVector();
k = k(1);
porosity = mat.getMaterial(1).PorousRock.getPorosity();
dyn = mat.getFluid().getDynViscosity;
cf = mat.getFluid().getFluidCompressibility; % [kPa^-1] Fluid compressibility
E = mat.getMaterial(1).ConstLaw.E;
nu = mat.getMaterial(1).ConstLaw.nu;
biot = mat.getMaterial(1).PorousRock.getBiotCoefficient();

% compute material parameters for analytical formulas
lambda =(E*nu)/((1+nu)*(1-2*nu)); %[kPa] first lamè constant
G = E/(2*(1+nu)); %[kPa] second lamè constant
M = (porosity*cf)^-1; %Biot Modulus, assuming cbr=0
Ku = lambda + 2*(G/3) + biot^2*M;
B = biot*M/Ku;
nu = lambda/(2*(lambda+G));
nuU = (3*nu+biot*B*(1-2*nu))/(3-biot*B*(1-2*nu)); %undrained Poisson Coefficient
c = 2*k*G*(1-nu)*(nuU-nu)/(dyn*biot^2*(1-nuU)*(1-2*nu)^2);
gamma = B*(1+nuU)/(3*(1-nuU));

% nz = 50; %number of calculation points along z-axis
% z = linspace(L_min,L_max,nz);

zu = mesh.coordinates(:,3);

if isfield(grid,'faces')
  zp = mesh.cellCentroid(:,3);
else
  zp = zu;
end
time = solv.domains.outstate.timeList;

zu = reshape(zu,1,[]);
zp = reshape(zp,1,[]);
time = reshape(time,1,[]);

[~,u0] = iniSol(zu,zp,M,F,Ku,biot,G);
[p,u] = TerzaghiSol(u0,zu,zp,time,nm,L,c,F,biot,gamma,G,nu);

% fileOut = fullfile(outDir,"Terzaghi_Analytical.mat");
% save(fileOut,'p','u','z','time')
% fprintf('Done computing Terzaghi analytical solution.\n');
end

function [p,u] = TerzaghiSol(u0,zu,zp,t,nm,L,c,pL,biot,gamma,G,nu)

nzu = length(zu);
nzp = length(zp);
m = 0:nm;
nt = length(t);

cellsP = arrayfun(@(m) bsxfun(@(zp,t) (1/(2*m+1))*exp((-(2*m+1)^2*pi^2*c*t)/(4*L^2))*sin(((2*m+1)*pi*(L-zp))/(2*L)),zp',t),m,'UniformOutput',false);
cellsU = arrayfun(@(m) bsxfun(@(zu,t) (1/(2*m+1)^2)*exp((-(2*m+1)^2*pi^2*c*t)/(4*L^2))*cos(((2*m+1)*pi*(L-zu))/(2*L)),zu',t),m,'UniformOutput',false);

matP = cell2mat(cellsP);
matU = cell2mat(cellsU);

seriesP = sum(reshape(matP,nzp,nt,[]),3);
seriesU = sum(reshape(matU, nzu,nt,[]),3);

p = (4*gamma*pL/pi)*seriesP; %[kPa]
u = repmat(u0',1,nt) + ((1-2*nu)*biot*gamma*pL/2/G/(1-nu))*(repmat(zu',1,nt) - (8*L/pi^2)*seriesU); 

end

function [p0,u0] = iniSol(zu,zp,M,pL,Ku,biot,G)
nzp = length(zp);
p0 = zeros(nzp,1);
p0(1:end) = (biot*M*pL)/(Ku+4*G/3);
u0 = arrayfun(@(zu) 1/(Ku+4*G/3)*pL*(zu),zu);
end



