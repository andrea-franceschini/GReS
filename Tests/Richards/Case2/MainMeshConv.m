close all;
% clear;
output_dir = 'Outputs';
input_dir = 'Inputs';
figures_dir = fullfile(output_dir,"Images");

%% ------------------------------ Set up the Domain -----------------------
% Set the simulation parameters for the non-linear solver.
simParam = SimulationParameters(fullfile(input_dir,'simparamMeshConv.xml'));

% Create an object of the Materials class and read the materials file
mat = Materials(fullfile(input_dir,'materials.xml'));

result = struct('pressure',[],'saturation',[],'height',[]);
meshList = [10 20 40 80 160];
for i=1:length(meshList)
    % Create the Mesh object
    topology = Mesh();

    % Set the input file name
    fileName = fullfile(input_dir,'Mesh',strcat('Column',int2str(meshList(i)),'.msh'));

    % Import mesh data into the Mesh object
    topology.importMesh(fileName);

    % Create an object of the "Elements" class and process the element properties
    elems = Elements(topology,2);

    % Create an object of the "Faces" class and process the face properties
    faces = Faces(topology);

    % Wrap Mesh, Elements and Faces objects in a structure
    grid = struct('topology',topology,'cells',elems,'faces',faces);

    % Creating boundaries conditions.
    bound = Boundaries(fullfile(input_dir,'boundaries.xml'),grid);

    % ------------------ Set up and Calling the Solver -----------------------
    % Create and set the print utility
    printUtils = OutState(topology,fullfile(input_dir,'outputMeshConv.xml'));

    % Create object handling construction of Jacobian and rhs of the model
    domain = Discretizer('Grid',grid,...
                         'Materials',mat,...
                         'Boundaries',bound,...
                         'OutState',printUtils);

    domain.addPhysicsSolver(fullfile(input_dir,'solver.xml'));

    % set initial conditions directly modifying the state object
    domain.state.data.pressure(:) = -9.8066e4;

    % The modular structure of the discretizer class allow the user to easily
    % customize the solution scheme.
    % Here, a built-in fully implict solution scheme is adopted with class
    % FCSolver. This could be simply be replaced by a user defined function
    % Solver = FCSolver(domain);
    Solver = GeneralSolver(simParam,domain);

    % Solve the problem
    Solver.NonLinearLoop();

    % Finalize the print utility
    printUtils.finalize()

    % Retriving the analisys information
    pressure = printUtils.results(1).pressure;
    saturation = printUtils.results(1).saturation;

    %Getting pressure and saturation solution for specified time from MatFILE
    numb = 0.;
    tol = 0.01;
    nodesP = find(abs(elems.mesh.cellCentroid(:,1)-numb) < tol & abs(elems.mesh.cellCentroid(:,2)-numb) < tol);

    result(i).height = elems.mesh.cellCentroid(nodesP,3);
    result(i).pressure = pressure(nodesP);
    result(i).saturation = saturation(nodesP);

    clearvars -except input_dir figures_dir model simParam mat result meshList i
end

%% --------------------- Post Processing the Results ----------------------
image_dir = fullfile(pwd,figures_dir);
if ~isfolder(image_dir)
    mkdir(image_dir)
end

% Ajusting the legend
tstr = strcat(int2str(meshList(:)),' Cells');

% Ajusting the pressure
weight = mat.getFluid().getSpecificWeight();

%Plotting pressure head
figure('Position', [100, 100, 700, 700])
hold on
for i=1:length(result)
    plot(result(i).pressure/weight,result(i).height,'.-', 'LineWidth', 2, 'MarkerSize', 14);
end
xlabel('Pressure (m)')
ylabel('Height (m)')
legend(tstr, 'Location', 'northwest')
set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
% export figure with quality
stmp = fullfile(image_dir,'mesh_pressure.png');
exportgraphics(gcf,stmp,'Resolution',400)

%Plotting saturation
figure('Position', [100, 100, 700, 700])
hold on
for i=1:length(result)
    plot(result(i).saturation,result(i).height,'.-', 'LineWidth', 2, 'MarkerSize', 14);
end
xlabel('Saturation')
ylabel('Height (m)')
legend(tstr, 'Location', 'northwest')
str = strcat('t = ',tstr);
set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
% export figure with quality
stmp = fullfile(image_dir,'mesh_saturation.png');
exportgraphics(gcf,stmp,'Resolution',400)
