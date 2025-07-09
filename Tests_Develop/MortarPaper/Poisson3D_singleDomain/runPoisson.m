% Import the mesh data into the Mesh object
topology.importMesh(fullfile(strcat(fname,'.vtk')));

% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,nG);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);
%

% Degree of freedom manager 
%fname = 'dof.dat';
dofmanager = DoFManager(topology,model);

% Create object handling construction of Jacobian and rhs of the model
linSyst = Discretizer(model,simParam,dofmanager,grid,mat);

linSyst.getSolver('Poisson').setAnalSolution(anal,f,gradx,grady,gradz);

% Build a structure storing variable fields at each time step
linSyst.setState();

% Create and set the print utility
printUtils = OutState(model,topology,'outTime.dat','folderName',strcat('out',num2str(i)));


% write files for bcs
nodes = unique(topology.surfaces);
c = topology.coordinates(nodes,:);
vals = arrayfun(@(i) anal(c(i,:)),1:numel(nodes));
vals = reshape(vals,[],1);
writeBCfiles('bc','NodeBC','Dir','Poisson','manufactured_bc',0,0,nodes,vals);

% Create an object of the "Boundaries" class 
bound = Boundaries("bc.dat",model,grid);

% Print model initial state
printState(printUtils,linSyst);

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme. 
% Here, a built-in fully implict solution scheme is adopted with class
% FCSolver. This could be simply be replaced by a user defined function
Solver = FCSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,linSyst);
%
% Solve the problem
[simState] = Solver.NonLinearLoop();
%
% Finalize the print utility
printUtils.finalize()
