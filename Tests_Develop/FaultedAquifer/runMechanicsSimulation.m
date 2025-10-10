function runMechanicsSimulation(mesh,pressures)

fprintf('Faulted aquifer model - Mechanical simulation \n')
fprintf('___________________\n\n')


% Set physical models 
model = ModelType("Poromechanics_FEM");

% Set parameters of the simulation
fileName = "simParam_mech.dat";
simParam = SimulationParameters(fileName,model);


% Create an object of the Materials class and read the materials file
fileName = 'Materials/materialsList.dat';
mat = Materials(model,fileName);


% Create an object of the "Elements" class and process the element properties
gaussOrder = 2;
elems = Elements(mesh,gaussOrder);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, mesh);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',mesh,'cells',elems,'faces',faces);
%
dofmanager = DoFManager(mesh, model);

% Create and set the print utility
printUtils = OutState(model,mesh,'outTime.dat','folderName','Faulted_aquifer_mech','flagMatFile',false);

% Create an object of the "Boundaries" class 
% write BC files
setMechBCfiles(mesh,pressures);
fileName = ["BCs/pressures.dat", "BCs/fix_X.dat", "BCs/fix_Y.dat", "BCs/fix_Z.dat"];
bound = Boundaries(fileName,model,grid);


% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('ModelType',model,...
                     'SimulationParameters',simParam,...
                     'DoFManager',dofmanager,...
                     'Boundaries',bound,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);

interfFile = 'Domains/interface.xml';

% set verbosity 
[interfaces,domain] = Mortar.buildInterfaces(interfFile,domain);

%setInitialTraction(interfaces{1});

interfaces{1}.solvers(2).simparams.setVerbosity(2);

solver = MultidomainFCSolver(domain,interfaces);
%maxActiveSetIters = 10;
% interfaces{1}.
%solver = ActiveSetContactSolver(domain,interfaces,maxActiveSetIters);

%solver.simParameters.setBackstepSkipFlag(1);

solver.NonLinearLoop();
solver.finalizeOutput();
end




function setMechBCfiles(mesh,press)
mkdir BCs
% custom BCs

% fix displacements in outer boundaries
writeBCfiles('BCs/fix_Z','SurfBC','Dir',{'Poromechanics','z'},'fix_Z',0,0,mesh,3);
writeBCfiles('BCs/fix_Y','SurfBC','Dir',{'Poromechanics','y'},'fix_Y',0,0,mesh,[4 6]);
writeBCfiles('BCs/fix_X','SurfBC','Dir',{'Poromechanics','x'},'fix_X',0,0,mesh,[5 7]);

% write pressure bcs

times = [0,3,7,10];
press = 1e-3*press(:,times+1);

list = (1:mesh.nCells)';

fName = 'BCs/pressures';
if ~isfolder(fName)
  mkdir(fName);
end
listName = strcat(fName,'/list');
fList = fopen(listName,'w');


fprintf(fList,'%i         %% Number of fixed entities \n',length(list));
fprintf(fList,'%i \n',list);

for i = 1:length(times)
  vals = press(:,i);
  t_name = strcat(fName,'/time',num2str(i-1),'.dat');
  ft = fopen(t_name,'w');
  fprintf(ft,'%%Time %2.4f \n',times(i));
  fprintf(ft,'%1.6e \n',vals);
end

end


function setInitialTraction(interface)

K0 = 1-sin(deg2rad(30)); % horizontal factor

gamma_s = 0.2; %specific weight of soil

sigma_v = gamma_s*interface.mesh.msh(2).surfaceCentroid(:,3);

sigma_glob = [-K0*sigma_v -K0*sigma_v -sigma_v];

for i = 1:interface.mesh.msh(2).nSurfaces
  s = diag(sigma_glob(i,:));
  R = interface.contactHelper.getRotationMatrix(i);
  n = R(:,1);
  sloc = R'*(s*n); % Rt*(sigma*n)
  interface.multipliers.curr(dofId(i,3)) = sloc;
  interface.iniMultipliers(dofId(i,3)) = sloc;
  interface.multipliers.prev(dofId(i,3)) = sloc;
end

end