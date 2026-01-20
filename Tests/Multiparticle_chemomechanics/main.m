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

%% BUILD MODEL

simparams = SimulationParameters(fullfile(input_dir, ...
    'simparam.xml'));
% Initialize the mortar utilities
[domains, interfaces] = buildModel(fullfile(input_dir, ...
    'multiparticle_chemomechanics.xml'));


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
