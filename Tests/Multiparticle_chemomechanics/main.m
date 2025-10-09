close all
clear

% Set some model parameters not included in the parameters file
c_max = 2.29e4;                     % maximum concentration [mol.m-3]
D = 7.08e-15;                        % diffusion coefficient [m2.s-1]
F = 96485;                          % Faraday constant [C.mol-1]
Rp = 5e-6;                          % particle radius [m]
Crate = 1;                          % C-rate for charge
i_n = F*Crate/3600 * c_max * Rp/3;  % Normal current density [A.m-2]
I = i_n*Rp / (D*c_max*F);           % Nondimensional normal current density

%% PREPROCESSING TO WRITE BCS FILE PROGRAMMATICALLY

% Load the mesh for setting up boundary conditions
meshFile = fullfile('Mesh','two_spheres.vtk');
mesh = Mesh();
mesh.importMesh(meshFile);

% get dirichlet node index (X = MIN, Y = MAX, Z = MAX)
tol = 1e-3;

[cx,cy,cz] = deal(mesh.coordinates(:,1),...
                  mesh.coordinates(:,2),...
                  mesh.coordinates(:,3));
ContactPoint = find(all(abs([cx, cy cz]) < tol,2));
mesh.surfaceTag(ContactPoint) = NaN; % Remove ContactPoint from surfaces
writeBCfiles('BCs/LeftSphere_gal', 'SurfBC', 'Neu', 'SinglePhaseFlow', ...
    'Left_1C', 0, I, 3);
writeBCfiles('BCs/RightSphere_gal', 'SurfBC', 'Neu', 'SinglePhaseFlow', ...
    'Right_1C', 0, I, 4);
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

%% RUN MODEL  
% A different solver is needed for models with non conforming domains
solver = MultidomainFCSolver(domains,interfaces);
solver.NonLinearLoop();
solver.finalizeOutput();
