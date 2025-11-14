close all
clear
tic;

% Set some model parameters not included in the parameters file
material_switch = 0; % 0 for [Zhang,2007], 1 for Si
% Make sure to also switch Omega_d in Poromechanics.m
if material_switch == 0
    c_max = 2.29e4;                     % maximum concentration [mol.m-3]
    D = 7.08e-15;                       % diffusion coefficient [m2.s-1]
    Rp = 5e-6;                          % particle radius [m]
elseif material_switch == 1
    c_max = 2.95e5;                     % maximum concentration [mol.m-3]
    D = 1.0e-16;                        % diffusion coefficient [m2.s-1]
    Rp = 500e-9;                        % particle radius [m]
else
    error(['Please provide a valid material_switch (0 for [Zhang,2007]', ...
        'and 1 for Si)']);
end

F = 96485;                          % Faraday constant [C.mol-1]
Crate = 1;                          % C-rate for charge
i_n = F*Crate/3600 * c_max * Rp/3;  % Normal current density [A.m-2]
I = i_n*Rp / (D*c_max*F);           % Nondimensional normal current density

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
    'Left_1C', 0, I, mesh, 3);
writeBCfiles('BCs/RightSphere_gal', 'SurfBC', 'Neu', 'SinglePhaseFlow', ...
    'Right_1C', 0, I, mesh, 4);

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

%% RUN MODEL  
% A different solver is needed for models with non conforming domains
solver = MultidomainFCSolver(domains,interfaces);
solver.NonLinearLoop();
solver.finalizeOutput();
toc;
