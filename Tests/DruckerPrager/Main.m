close all;
clear;
clc;

scriptFullPath = mfilename('fullpath');
scriptDir = fileparts(scriptFullPath);
cd(scriptDir);

% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("Poromechanics_FEM");

% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = "simParam.dat";
simParam = SimulationParameters(fileName,model);

% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();

% Set the input file name
fileName = 'Mesh/cilinder.msh';

% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);

%----------------------------- MATERIALS -----------------------------
% Set the input file name
fileName = 'Materials/MaterialsList.dat';

% Create an object of the Materials class and read the materials file
mat = Materials(model,fileName);

%------------------------------ ELEMENTS -----------------------------
GaussPts = Gauss(12,2,3);

% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,GaussPts);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);

% Wrap Mesh, Elements and Faces objects in a structure
grid1 = struct('topology',topology,'cells',elems,'faces',faces);

% Create a DoFManager object
dofmanager = DoFManager(topology, model);

% Create object handling construction of Jacobian and rhs of the model
linSyst = Discretizer(model,simParam,dofmanager,grid1,mat,GaussPts);

% Build a structure storing variable fields at each time step
state = linSyst.setState();

% Create and set the print utility
printUtils = OutState(model,topology,'outTime.dat','folderName','vtkOutput');

%------------------------ BOUNDARY CONDITIONS ------------------------
% Set the input file
fileName = ["BCs/bottom_fix_cilinder2.dat", "BCs/load_z_cilinder2.dat"];    

% Create an object of the "Boundaries" class and read the boundary conditions
bound = Boundaries(fileName,model,grid1);

% Set initial stress state
nu = mat.db(1).ConstLaw.nu;
z0 = -60;   % average depth
gamma = 120; % N/m^3
M1 = 0.5;
M2 = 0.5;
cells = dofmanager.getFieldCells("Poromechanics");
for el = cells'
  switch topology.cellVTKType(el)
    case 10
      z1 = z0 + elems.cellCentroid(el,3);
      sigmaz = gamma*z0;
      sigmax = M1*sigmaz;
      sigmay = M2*sigmaz;
      state.iniStress(el,1) = sigmax;
      state.iniStress(el,2) = sigmay;
      state.iniStress(el,3) = sigmaz;
  end
end
state.conv.stress = state.iniStress;

% ---------------------------- SOLUTION -------------------------------
% Create the object handling the (nonlinear) solution of the problem
solver = FCSolver(model,simParam,dofmanager,grid1,mat,bound,printUtils,state,linSyst,GaussPts);
[simState, endState] = solver.NonLinearLoop();

% Finalize the print utility
printUtils.finalize()

%% POST PROCESSING

% Initialize and plot results
image_dir = strcat(pwd,'/Images');
if isfolder(image_dir)
    rmdir(image_dir,"s")
    mkdir Images
else
    mkdir Images
end

% Find cells near axis
tol = 0.005;
elem_axis_x = find(abs(elems.cellCentroid(:,1))<tol);
elem_axis_y = find(abs(elems.cellCentroid(:,2))<tol);
elem_axis = intersect(elem_axis_x,elem_axis_y);

% Sort elements along z-axis
elem_z_coords = elems.cellCentroid(elem_axis, 3);
[elem_z_coords_sorted, sort_idx] = sort(elem_z_coords);
elem_axis_sorted = elem_axis(sort_idx);

load("resultsDP.mat");

% Extract stress history for elements on axis
n_elements = length(elem_axis_sorted);
n_timesteps = 20; 

numerical_stress_history = cell(n_elements, 1);
elem_coords_for_analytical = elem_z_coords_sorted;

% Extract stress data for each element and each timestep
for i = 1:n_elements
    elem_idx = elem_axis_sorted(i);
    numerical_stress_history{i} = cell(n_timesteps, 1);
    
    for j = 1:n_timesteps
        stress_data = stress_story{elem_idx}{j};
        
        if isstruct(stress_data)
            numerical_stress_history{i}{j} = stress_data.stress(:);
        elseif length(stress_data) >= 3
            if length(stress_data) == 3
                numerical_stress_history{i}{j} = [stress_data(:); 0; 0; 0];
            else
                numerical_stress_history{i}{j} = stress_data(:);
            end
        else
            error('stress_story format not recognized for element %d, step %d', elem_idx, j);
        end
    end
end

% Analytical solution
fprintf('Computing analytical solution...\n');
DruckerPrager_analytical(mat, z0, gamma, M1, M2, ...
                                numerical_stress_history, elem_coords_for_analytical);

% Load analytical results
load("DruckerPrager_Analytical.mat");
load('StrainResults.mat');

elem_idx_num = 2;
elem_idx_anal = 1;

% --- Analytical ---
eps_z_analytical = squeeze(epsilon_total(elem_idx_anal, 3, :))';

% --- Numerical ---
n_steps = length(strain_story{elem_idx_num});
eps_z_numerical = zeros(1, n_steps);

% Time vector
t = simParam.dtIni:simParam.dtMax:simParam.dtIni + (length(strain_story{elem_axis_sorted(1)})-1)*simParam.dtMax;

% Choose step to analyze (>=2)
n_steps = length(eps_z_numerical);
step_k = 20;
if step_k < 2
    step_k = 2;
end

% Calculate increments (Δx_k = x_k - x_{k-1})
d_eps_zz_num = [0, diff(eps_z_numerical)];
d_eps_zz_an = [0, diff(eps_z_analytical)];

% Difference between increments (Δnum − Δanal) at each step
diff_incr_zz = d_eps_zz_num - d_eps_zz_an;

err_start_zz = abs(eps_z_numerical(step_k-1) - eps_z_analytical(step_k-1));
err_end_zz = abs(eps_z_numerical(step_k) - eps_z_analytical(step_k));
mismatch_zz = diff_incr_zz(step_k);

% Plot differences
fig_inc = figure('Name','Difference between analytical and numerical solution', 'Position',[300 300 1000 600]);
plot([step_k-1, step_k], [err_start_zz, err_end_zz], 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Timestep k'); 
ylabel('\epsilon_{zz,num} - \epsilon_{zz,anal}');
title('Difference between \epsilon_{zz} of analytical and numerical solution');
ylim_curr = ylim;
ylim([min(0, ylim_curr(1)), max(0, ylim_curr(2))]);
exportgraphics(fig_inc, fullfile(image_dir, 'Difference between analytical and numerical solution.png'), 'Resolution', 300);

fprintf('[STEP: %d]  eps_zz:  Δnum−Δanal = %.6e | initial difference = %.6e | final difference = %.6e\n', ...
        step_k, mismatch_zz, err_start_zz, err_end_zz);