function runContactMechanicsSimulation(fileName,mesh,pressures)

fprintf('Faulted aquifer model - Contact simulation \n')
fprintf('___________________________________________\n\n')

simParam = SimulationParameters(fileName);

mat = Materials(fileName);

% Create an object of the "Elements" class and process the element properties
gaussOrder = 2;
elems = Elements(mesh,gaussOrder);

% Create an object of the "Faces" class and process the face properties
faces = Faces(mesh);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',mesh,'cells',elems,'faces',faces);

% write BC values
setMechBCfiles(mesh,pressures);

bound = Boundaries(fileName,grid);

printUtils = OutState(mesh,fileName);


% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Boundaries',bound,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);

domain.addPhysicsSolver(fileName);


% set verbosity
interfaces = buildInterfaces(fileName,domain);

initTrac = setInitialTraction(interfaces{1});

if contains(fileName,"stick")
  % interfaces{1}.state.multipliers = initTrac;
  % interfaces{1}.state.iniMultipliers = initTrac;
  % interfaces{1}.stateOld.multipliers = initTrac;
  % interfaces{1}.stateOld.iniMultipliers = initTrac;
else
  interfaces{1}.state.iniTraction = initTrac;
  interfaces{1}.state.traction = initTrac;
  interfaces{1}.stateOld.iniTraction = initTrac;
  interfaces{1}.stateOld.traction = initTrac;
end

% interfaces{1}.contactHelper.forceStickBoundary = true;
% interfaces{1}.contactHelper.tol.normalGap = 1e-5;
% interfaces{1}.contactHelper.tol.normalTrac = -1e-3;
% interfaces{1}.contactHelper.tol.slidingCheck = 5e-3;
% interfaces{1}.contactHelper.tol.minLimitTraction = 1e-4;  % below this value, the limit traction is set to 0
% interfaces{1}.contactHelper.tol.areaTol = 2e-3;
%interfaces{1}.contactHelper.resetActiveSet = true;

if ~contains(fileName,"stick")
  maxActiveSetIters = 12;
  solver = ActiveSetContactSolver(simParam,domain,interfaces,maxActiveSetIters);
else
  solver = MultidomainFCSolver(simParam,domain,interfaces);
end


solver.NonLinearLoop();
solver.finalizeOutput();

end




function setMechBCfiles(mesh,press)

% write pressure bcs

times = [1,3,7,10];

% convert pressure from kPa to MPa
press = 1e-3*press(:,times);

list = (1:mesh.nCells)';

fName = 'InputMech/pressureBC';
if ~isfolder(fName)
  mkdir(fName);
end

listName = strcat(fName,'/list');
fList = fopen(listName,'w');

% print entitiy file
fprintf(fList,'%i         %% Number of fixed entities \n',length(list));
fprintf(fList,'%i \n',list);

% print value file
for i = 1:length(times)
  vals = press(:,i);
  t_name = strcat(fName,'/time',num2str(i-1),'.dat');
  ft = fopen(t_name,'w');
  fprintf(ft,'%%Time %2.4f \n',times(i));
  fprintf(ft,'%1.6e \n',vals);
end

end


function iniTrac = setInitialTraction(interface)

K0 = 1-sin(deg2rad(30)); % horizontal factor

gamma_s = 0.021; %specific weight of soil

mshSlave = interface.getMesh(MortarSide.slave);

depth = abs(max(mshSlave.surfaceCentroid(:,3)) - mshSlave.surfaceCentroid(:,3));

coes = 0.12;
sigma_v = coes+gamma_s*depth;

sigma_glob = [-K0*sigma_v -K0*sigma_v -sigma_v];

iniTrac = zeros(interface.nMult,1);

for i = 1:mshSlave.nSurfaces
  s = diag(sigma_glob(i,:));
  R = interface.interfMesh.getRotationMatrix(i);
  n = R(:,1);
  sloc = R'*(s*n); % Rt*(sigma*n)
  dofs = DoFManager.dofExpand(i,3);
  iniTrac(dofs) = sloc;
end

end
