close all;
clear;
tic;
% Set some model parameters not included in the parameters file
c_max = 2.29e4;                     % maximum concentration [mol.m-3]
D = 7.08e-15;                        % diffusion coefficient [m2.s-1]
F = 96485;                          % Faraday constant [C.mol-1]
Rp = 5e-6;                          % particle radius [m]
Crate = 1;                          % C-rate for charge
i_n = F*Crate/3600 * c_max * Rp/3;  % Normal current density [A.m-2]
I = i_n*Rp / (D*c_max*F);           % Nondimensional normal current density

% Set physical models 
model = ModelType(["SinglePhaseFlow_FEM", "Poromechanics_FEM"]);

% Set parameters of the simulation
fileName = "simParam.dat";
simParam = SimulationParameters(fileName,model);

% Create the Mesh object
topology = Mesh();

% Set the mesh input file name
fileName = 'Mesh/3d_sphere.msh';
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);

% Create an object of the Materials class and read the materials file
fileName = 'materialsList_Zhang2007.dat';
mat = Materials(model,fileName);

% Create an object of the "Elements" class and process the element properties
gaussOrder = 1;
elems = Elements(topology,gaussOrder);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);

% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);

% Degree of freedom manager
dofmanager = DoFManager(topology,model);

% Create and set the print utility
printUtils = OutState(model, topology, 'outTime.dat', 'folderName', ...
    'Output_chemomech_tetra', 'flagMatFile', true);

% Write BC files programmatically with function utility
% Fixing the particle center (ux=uy=uz=0)
r = sqrt(sum(topology.coordinates.^2, 2)); % distance of each node from origin
r0 = 0.12 * max(r); % threshold constraint (e.g., 5% of sphere radius)
constrainIndices = find(r < r0); % all nodes within that small sphere
writeBCfiles('BCs/chemomech_u_0', 'NodeBC', 'Dir', {'Poromechanics', ...
    'x', 'y', 'z'}, 'Fixed center point', 0, 0, constrainIndices);
% Outer boundary condition is natural (normal stress = 0)

% Potentiostatic boundary condition (constant c)
writeBCfiles('BCs/chemomech_cmax', 'SurfBC', 'Dir', 'SinglePhaseFlow', ...
    'c_outer_bc', 0, 1, topology, 2); % 1 for nondimensional bc(=c_max)

% Galvanostatic boundary condition (constant I)
writeBCfiles('BCs/chemomech_gal', 'SurfBC', 'Neu', 'SinglePhaseFlow', ...
    'c_outer_bc', 0, I, topology, 2);

% Collect BC input file in a list
fileName = ["BCs/chemomech_gal.dat", "BCs/chemomech_u_0.dat"];

% Create an object of the "Boundaries" class
bound = Boundaries(fileName,model,grid);

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('ModelType',model,...
                     'SimulationParameters',simParam,...
                     'DoFManager',dofmanager,...
                     'Boundaries',bound,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);

% Apply initial conditions
c_in = 6195; % initial concentration value
c_in_d = c_in / c_max; % nondimensional value of initial concentration
domain.state.data.pressure = domain.state.data.pressure + c_in_d;
% Initial displacements are zero by default - no need to change

% Print model initial state
printState(domain);

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme. 
% Here, a built-in fully implict solution scheme is adopted with class
% FCSolver. This could be simply be replaced by a user defined function
Solver = FCSolver(domain);
%
% Solve the problem
[simState] = Solver.NonLinearLoop();
%
% Finalize the print utility
domain.outstate.finalize()
toc;