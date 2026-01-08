% close all
clear
tic;

% Add path to helper functions for stress/strain computation
scriptFullPath = mfilename('fullpath');
scriptDir = fileparts(scriptFullPath);
chemoMechDir = fullfile(fileparts(scriptDir), 'Chemo_mechanics');
if exist(chemoMechDir, 'dir')
    addpath(chemoMechDir);
end

% Set some model parameters not included in the parameters file
material_switch = 0                                               ; % 0 for [Zhang,2007], 1 for Si
% Make sure to also switch material files in domain.xml
if material_switch == 0
    params.c_max = 2.29e4;                     % maximum concentration [mol.m-3]
    params.D = 7.08e-15;                       % diffusion coefficient [m2.s-1]
    params.Rp = 5e-6;                          % particle radius [m]
    params.E = 10e9;            % elastic modulus [N.m-2] (for [Zhang,2007])
    params.nu = 0.3;            % Poisson's ratio
    params.Omega = 3.497e-6;    % partial molar volume [mol.m-3]
elseif material_switch == 1
    params.c_max = 2.95e5;                     % maximum concentration [mol.m-3]
    params.D = 1.0e-16;                        % diffusion coefficient [m2.s-1]
    params.Rp = 500e-9;                          % particle radius [m]
    params.E = 160e9;            % elastic modulus [N.m-2] (for Si - adjust as needed)
    params.nu = 0.22;           % Poisson's ratio (for Si)
    params.Omega = 8.89e-6;    % partial molar volume [mol.m-3] (adjust as needed)
else
    error(['Please provide a valid material_switch (0 for [Zhang,2007]', ...
        'and 1 for Si)']);
end

params.F = 96485;                          % Faraday constant [C.mol-1]
params.Crate = 1;                          % C-rate for charge
params.i_n = 2; % params.F*params.Crate/3600 * params.c_max * params.Rp/3;  % Normal current density [A.m-2]
params.I = params.i_n*params.Rp / (params.D*params.c_max*params.F);           % Nondimensional normal current density

% Nondimensional material parameters
params.E_d = params.E / params.E;  % = 1 (nondimensional)
params.G_d = params.E_d / (2*(1 + params.nu));
params.lambda_d = params.nu*params.E_d/((1 + params.nu)*(1 - 2*params.nu));
params.Omega_d = params.Omega * params.c_max;
params.beta_d = params.Omega_d * (3*params.lambda_d + 2*params.G_d) / 3;
params.c_in = 0;                  % initial concentration value
params.c_in_d = params.c_in / params.c_max; % nondimensional value of initial concentration

%% PREPROCESSING TO WRITE BCS FILE PROGRAMMATICALLY

% Load the mesh for setting up boundary conditions
meshFile = fullfile('Mesh','two_spheres.vtk');
mesh = Mesh();
mesh.importMesh(meshFile);

% get dirichlet node index (X = MIN, Y = MAX, Z = MAX)
tol = 1e-2;

% Fix the left sphere
[cx,cy,cz] = deal( mesh.coordinates(:,1), ...
                   mesh.coordinates(:,2), ...
                   mesh.coordinates(:,3) ...
                   );

% Poromechanics boundary condition (fix sphere centers)
[XL, YL, ZL] = deal(-1, 0, 0);
LeftSphereCenter = find(all(abs([cx-XL, cy-YL cz-ZL]) < tol,2));
writeBCfiles('BCs/FixLeftSphere', 'NodeBC', 'Dir', {'Poromechanics', ...
    'x', 'y', 'z'}, 'FixLeftSphereCenter', 0, 0, LeftSphereCenter);

[XR, YR, ZR] = deal(1, 0, 0);
RightSphereCenter = find(all(abs([cx-XR, cy-YR cz-ZR]) < tol,2));
writeBCfiles('BCs/FixRightSphere', 'NodeBC', 'Dir', {'Poromechanics', ...
    'x', 'y', 'z'}, 'FixRightSphereCenter', 0, 0, LeftSphereCenter);

[XC, YC, ZC] = deal(0, 0, 0);
ContactPoint = find(all(abs([cx-XC cy-YC cz-ZC]) < tol, 2));
d2 = cx(:).^2 + cy(:).^2 + cz(:).^2;   % squared distances (Nx1)

% exclude the point with index pt
d2(ContactPoint) = Inf;

[~, ContactPoint2] = min(d2);
writeBCfiles('BCs/FixContactPoint', 'NodeBC', 'Dir', {'Poromechanics', ...
    'x', 'y', 'z'}, 'FixContactPoint', 0, 0, ContactPoint);

% SinglePhaseFlow boundary conditions (galvanostatic)
writeBCfiles('BCs/LeftSphere_gal', 'SurfBC', 'Neu', 'SinglePhaseFlow', ...
    'Left_1C', 0, params.I, mesh, 3);
writeBCfiles('BCs/RightSphere_gal', 'SurfBC', 'Neu', 'SinglePhaseFlow', ...
    'Right_1C', 0, params.I, mesh, 4);

% SinglePhaseFlow boundary conditions (c = c_max)
writeBCfiles('BCs/LeftSphere_cmax', 'SurfBC', 'Dir', 'SinglePhaseFlow', ...
    'Left_cmax', 0, 1, mesh, 3);
writeBCfiles('BCs/RightSphere_cmax', 'SurfBC', 'Dir', 'SinglePhaseFlow', ...
    'Right_cmax', 0, 1, mesh, 4);
clear mesh

%% BUILD MODEL
% For now, simulation parameters are shared by all domains
simParam = SimulationParameters('simParam.dat');

% build model using domains input file (a shortcut to programmatically
% initialize separate model objects for each domain)
domainFile = fullfile('Domains','domain.xml');
domains = buildModel(domainFile);

% Initialize the mortar utilities
interfFile = fullfile('Domains','interfaces.xml');
[interfaces,domains] = Mortar.buildInterfaces(interfFile,domains);

% % write Bc logic and prepare files to manually redefine the boundaries
% domains.bcs = Boundaries();

% Apply initial conditions to each domain
for i = 1:length(domains)
    % Set initial pressure (concentration)
    domains(i).state.data.pressure = domains(i).state.data.pressure + params.c_in_d;
    % Ensure that pOld is getting updated with the pressure vector
    domains(i).state.data.pOld = domains(i).state.data.pressure;
    % Initial displacements are zero by default - no need to change
end

% Print model initial state for each domain
for i = 1:length(domains)
    printState(domains(i));
end

%% RUN MODEL  
% A different solver is needed for models with non conforming domains
solver = MultidomainFCSolver(domains,interfaces);
solver.NonLinearLoop();
solver.finalizeOutput();

%% POST-PROCESS STRESSES AND STRAINS
% Extract p and u values from printUtils for each domain
% Note: In this case, there's a single domain, but the code handles multiple domains
% IMPORTANT: Ensure that OutState has matFile="1" in domain.xml to store results
for i = 1:length(domains)
    if ~isempty(domains(i).outstate) && ~isempty(domains(i).outstate.results)
        % Check if results are populated
        if isempty([domains(i).outstate.results.expTime])
            warning('Domain %d: No results found. Ensure OutState has matFile="1" in domain.xml', i);
            continue;
        end
        
        % Extract output data
        output_times = [domains(i).outstate.results.expTime]'; % timesteps x 1
        p = [domains(i).outstate.results.expPress]'; % timesteps x nNodes
        u = [domains(i).outstate.results.expDispl]'; % timesteps x (3*nNodes)
        
        % Get mesh and elements from domain
        topology = domains(i).grid.topology;
        elems = domains(i).grid.cells;
        
        % Compute strains and stresses
        [strain, stress] = computeStrainsAndStresses_improved(output_times, ...
            p, u, topology, params, elems);
        
        % Get nodal strains and stresses
        strain_nodal = celltonodeStress(strain, topology);
        stress_nodal = celltonodeStress(stress, topology);
        
        fprintf('Domain %d: Computed strains and stresses for %d timesteps\n', ...
            i, length(output_times));
    else
        warning('Domain %d: outstate or results not available', i);
    end
end

run("plotting.m");

toc;
