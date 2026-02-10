% clear
close all;
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

%% Model parameters
% Set some model parameters not included in the parameters file
material_switch = 0; % 0 for [Zhang,2007], 1 for Si
% Make sure to also switch Materials file name in
% Inputs/multiparticle_chemomechanics.xml
if material_switch == 0
    params.c_max = 2.29e4;      % maximum concentration [mol.m-3]
    params.D = 7.08e-15;        % diffusion coefficient [m2.s-1]
    params.Rp = 5e-6;           % particle radius [m]
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

% get dirichlet node index (X = MIN, Y = MAX, Z = MAX)
tol = 1e-2;

%% BUILD MODEL

simparams = SimulationParameters(fullfile(input_dir, ...
    'simparam.xml'));
% Initialize the mortar utilities
[domains, interfaces] = buildModel(fullfile(input_dir, ...
    'multiparticle_chemomechanics.xml'));

% Extract mesh coordinates
[cx, cy, cz] = deal( domains.grid.topology.coordinates(:,1), ...
                   domains.grid.topology.coordinates(:,2), ...
                   domains.grid.topology.coordinates(:,3) ...
                   );

%% RUN MODEL  
% A different solver is needed for models with non conforming domains
solver = MultidomainFCSolver(simparams, domains, interfaces);
solver.NonLinearLoop();
solver.finalizeOutput();

%% POST-PROCESS STRESSES AND STRAINS
% Extract p and u values from printUtils for each domain
% Note: In this case, there's a single domain, but the code handles multiple domains
% IMPORTANT: Ensure that OutState has matFile="1" in domain.xml to store results
for i = 1:length(domains)
    if ~isempty(domains(i).outstate) && ~isempty(domains(i).outstate.results)
        % Check if results are populated
        if isempty([domains(i).outstate.results.time])
            warning('Domain %d: No results found. Ensure OutState has matFile="1" in domain.xml', i);
            continue;
        end
        
        % Extract output data
        output_times = [domains(i).outstate.results.time]'; % timesteps x 1
        p = [domains(i).outstate.results.pressure]'; % timesteps x nNodes
        u = [domains(i).outstate.results.displacements]'; % timesteps x (3*nNodes)
        
        % Get mesh and elements from domain
        topology = domains(i).grid.topology;
        elems = domains(i).grid.cells;
        
        % Compute strains and stresses
        [strain, stress] = computeStrainsAndStresses(output_times, ...
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
