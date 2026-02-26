close all;
clear;
tic;

%% Initialize input, output, directory paths
output_dir = 'Outputs';
input_dir = 'Inputs';
figures_dir = fullfile(output_dir,"Images");

% Add helper functions folder to path
addpath(genpath(fullfile(pwd, 'Helper functions/')));

% Extract the directory containing the script
scriptFullPath = mfilename('fullpath'); 
scriptDir = fileparts(scriptFullPath);

%% Simulation and model parameters
% Set parameters of the simulation
simParam = SimulationParameters(fullfile(scriptDir,input_dir,"simparam.xml"));

% Create an object of the Materials class and read the materials file
mat = Materials(fullfile(scriptDir,input_dir,"materials.xml"));

% Set some model parameters not included in the parameters file
params.E = 10e9;            % elastic modulus [N.m-2]
params.nu = 0.3;            % Poisson's ratio
params.D = 7.08e-15;        % diffusion coefficient [m2.s-1]
params.Omega = 3.497e-6;    % partial molar volume [mol.m-3]
params.c_max = 2.29e4;      % maximum concentration [mol.m-3]
params.F = 96485;           % Faraday constant [C.mol-1]
params.Rp = 5e-6;           % particle radius [m]

% % Normal current density - Crate specified [A.m-2]
% params.Crate = 1;
% params.i_n = params.F*params.Crate/3600 * params.c_max * params.Rp/3;

% Normal current density - i_n specified [A.m-2]
params.i_n = 2; % = 2 for [Zhang,2007] validation

% Nondimensional normal current density (I)
params.I = params.i_n*params.Rp / (params.D*params.c_max*params.F);

params.E_d = params.E / params.E;
params.G_d = params.E_d / (2*(1 + params.nu));
params.lambda_d = params.nu*params.E_d/((1 + params.nu)*(1 - 2*params.nu));
params.Omega_d = params.Omega * params.c_max;
params.beta_d = params.Omega_d * (3*params.lambda_d + 2*params.G_d) / 3;

%% Importing and setting up the mesh
% Importing mesh
topology = Mesh();
MeshFile = "two_overlapping_spheres.vtk";
topology.importMesh(fullfile(scriptDir, input_dir, "Mesh", MeshFile));

% Create an object of the "Elements" class and process the element properties
gaussOrder = 1;
elems = Elements(topology, gaussOrder);

% % Create an object of the "Faces" class and process the face properties
% faces = Faces(topology);

% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology', topology, 'cells', elems, 'faces', []);

%% Creating boundaries conditions.
% Extracting nodal coordinates
tol = 1e-3;
[cx,cy,cz] = deal( topology.coordinates(:,1), ...
                   topology.coordinates(:,2), ...
                   topology.coordinates(:,3) ...
                   );
r0 = 0.2*max(sqrt(cx.^2 + cy.^2 + cz.^2)); % threshold constraint

% Computing specific points from nodal coordinates
Point_Minus100 = find(all(abs([cx+1 cy-0 cz-0]) < tol, 2));
writematrix(Point_Minus100, "Inputs/Point_-100.dat");

% Remember to update current values and points in boundaries when changing
% i_n/C-rate values or when changing the mesh file.
bound = Boundaries(fullfile(scriptDir, input_dir, "boundaries.xml"), grid);

%% ------------------ Set up and Calling the Solver -----------------------
% Create and set the print utility for the solution
printUtils = OutState(topology, fullfile(scriptDir, input_dir, ...
    'output.xml'));

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Boundaries', bound,...
                     'OutState', printUtils,...
                     'Materials', mat,...
                     'Grid', grid);

domain.addPhysicsSolver('solver_FEM.xml');

% Apply initial conditions
params.c_in = 0; % 6195; % initial concentration value
params.c_in_d = params.c_in / params.c_max; % nondimensional value of initial concentration
domain.state.data.pressure = domain.state.data.pressure + params.c_in_d;

% Ensure that pOld is getting updated with the pressure vector
domain.state.data.pOld = domain.state.data.pressure;

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme. 
% Here, a built-in fully implict solution scheme is adopted with class
% FCSolver. This could be simply be replaced by a user defined function
Solver = FCSolver(simParam,domain);

% Solve the problem
[simState] = Solver.NonLinearLoop();

% Finalize the print utility
domain.outstate.finalize()

%% --------------------- Post Processing the Results ----------------------
if true
    image_dir = fullfile(pwd,figures_dir);
    if isfolder(image_dir)
    rmdir(image_dir,"s")
    mkdir(figures_dir)
    else
    mkdir(figures_dir)
    end
    
    output_times = [printUtils.results.time]'; % timesteps x 1
    p = [printUtils.results.pressure]'; % timesteps x nNodes
    u = [printUtils.results.displacements]'; % timesteps x (3*nNodes)
    
    [strain, stress] = computeStrainsAndStresses(output_times, ...
        p, u, topology, params, elems);
    
    % Get nodal strains and stresses
    strain_nodal = celltonodeStress(strain, topology);
    stress_nodal = celltonodeStress(stress, topology);
    
    % Plotting
    run("plotting.m");
end
toc;