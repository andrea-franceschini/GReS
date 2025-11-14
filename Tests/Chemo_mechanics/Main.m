close all;
clear;
tic;
% Set some model parameters not included in the parameters file
params.E = 10e9;            % elastic modulus [N.m-2]
params.nu = 0.3;            % Poisson's ratio
params.D = 7.08e-15;        % diffusion coefficient [m2.s-1]
params.Omega = 3.497e-6;    % partial molar volume [mol.m-3]
params.c_max = 2.29e4;      % maximum concentration [mol.m-3]
params.F = 96485;           % Faraday constant [C.mol-1]
params.Rp = 5e-6;           % particle radius [m]
params.Crate = 1;           % C-rate for charge
% Normal current density [A.m-2]
% params.i_n = params.F*params.Crate/3600 * params.c_max * params.Rp/3;
params.i_n = 2; % for [Zhang,2007] validation
% Nondimensional normal current density
params.I = params.i_n*params.Rp / (params.D*params.c_max*params.F);

params.E_d = params.E / params.E;
params.G_d = params.E_d / (2*(1 + params.nu));
params.lambda_d = params.nu*params.E_d/((1 + params.nu)*(1 - 2*params.nu));
params.Omega_d = params.Omega * params.c_max;
params.beta_d = params.Omega_d * (3*params.lambda_d + 2*params.G_d) / 3;

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
tol = 1e-1;
[cx,cy,cz] = deal( topology.coordinates(:,1), ...
                   topology.coordinates(:,2), ...
                   topology.coordinates(:,3) ...
                   );
r0 = 0.2*max(sqrt(cx.^2 + cy.^2 + cz.^2)); % threshold constraint

% Fix ux=uy=uz=0 at relevant points on the surface
% 1. Fix ux=0 at (0,0,1) and uy=uz=0 at (1,0,0)
    writeBCfiles_PointwiseConstraint(cx, cy, cz, tol);

% 2. Fix ux=uy=uz=0 at all points r<r0
    % writeBCfiles_FixCenterPoints(cx, cy, cz, r0);

% % Add the required Neumann boundary conditions for stress
biot = mat.getMaterial(1).PorousRock.getBiotCoefficient();
sigma_n_corrective = 0; % biot * 1; % as c_max=1 is imposed for SinglePhaseFlow
writeBCfile_sigma_nn('BCs/chemomech_sigma_n', 'Corrective_sigma_n', 0, ...
    -sigma_n_corrective, topology, 2);

% Potentiostatic boundary condition (constant c)
writeBCfiles('BCs/chemomech_cmax', 'SurfBC', 'Dir', 'SinglePhaseFlow', ...
    'c_outer_bc', 0, 1, topology, 2); % 1 for nondimensional bc(=c_max)

% Galvanostatic boundary condition (constant I)
writeBCfiles('BCs/chemomech_gal', 'SurfBC', 'Neu', 'SinglePhaseFlow', ...
    'c_outer_bc', 0, params.I, topology, 2);

% Collect BC input file in a list
% Case 1
fileName = ["BCs/chemomech_gal.dat", "BCs/chemomech_u_0.dat", ...
    "BCs/chemomech_u_x.dat", "BCs/chemomech_u_yz.dat", ...
    "BCs/chemomech_sigma_n.dat"];

% Case 2
% fileName = ["BCs/chemomech_gal.dat", "BCs/chemomech_u_0.dat", ...
%     "BCs/chemomech_sigma_n.dat"];
% fileName = ["BCs/chemomech_cmax.dat", "BCs/chemomech_u_0.dat"];

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
params.c_in = 0; % 6195; % initial concentration value
params.c_in_d = params.c_in / params.c_max; % nondimensional value of initial concentration
domain.state.data.pressure = domain.state.data.pressure + params.c_in_d;
% Ensure that pOld is getting updated with the pressure vector
domain.state.data.pOld = domain.state.data.pressure;

% Initial displacements are zero by default - no need to change

% Print model initial state
printState(domain);

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme. 
% Here, a built-in fully implict solution scheme is adopted with class
% FCSolver. This could be simply be replaced by a user defined function
Solver = FCSolver(domain);

% Solve the problem
[simState] = Solver.NonLinearLoop();

% Finalize the print utility
domain.outstate.finalize()

% Post-process strains and stresses from printUtils.results
output_times = [printUtils.results.expTime]'; % timesteps x 1
p = [printUtils.results.expPress]'; % timesteps x nNodes
u = [printUtils.results.expDispl]'; % timesteps x (3*nNodes)
[strain, stress] = computeStrainsAndStresses(output_times, p, u, ...
    topology, params);

% Get nodal strains and stresses
strain_nodal = celltonodeStress(strain, topology);
stress_nodal = celltonodeStress(stress, topology);

% Plotting
run("plotting.m");
toc;