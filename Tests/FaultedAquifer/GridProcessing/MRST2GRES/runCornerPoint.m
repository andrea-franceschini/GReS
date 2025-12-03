clear
clc

%% this script runs several test grid provided by MRST

% 1) Create a 20-by-20-by-5 faulted grid with a fault drop of 0.15 (m).
grdecl = simpleGrdecl([6, 5, 100], 0.15);

% Create the grid data structure.
G = processGRDECL(grdecl,'PreserveCpNodes',true,'Tolerance',0.01);

% Plot the resulting geometry.
figure (1), hg1 = plotGrid(G, 'FaceAlpha', 0.625);
view(3), grid on, axis tight


% rotate the view to see the fault structures from top

%% norne model
grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

% Save ACTNUM until later and then override the ACTNUM values to obtain a
% grid that contains all cells
actnum        = grdecl.ACTNUM;
grdecl.ACTNUM = ones(prod(grdecl.cartDims),1);
% G             = processGRDECL(grdecl,'checkgrid', false);

% Create the grid data structure.
G = processGRDECL(grdecl,'PreserveCpNodes',true,'Tolerance',0);

figure(1)

plotGrid(G,'FaceColor','none','EdgeAlpha',.1);
plotFaces(G, G.faces.tag>0, ...
'FaceColor','red','FaceAlpha',.2, ...
'EdgeColor','r','EdgeAlpha',.1);
axis off; view(0,90); zoom(1.7);

figure (2)
G = computeGeometry(G,'CpGeometry',true);
%G.nodes.coords(:,3) = G.nodes.coords(:,3)*10;
% get faces with 4 faces only
cell_tetra = find(diff(G.cells.facePos)==4);
plotGrid(G,'FaceColor','none','EdgeAlpha',.05);
plotGrid(G,cell_tetra,'FaceColor','y','EdgeAlpha',.5);
%%
figure (3)
id = G.cells.cpgeometry.extent(:,3) < 1e-3;
% get faces with 4 faces only
%cell_tetra = find(diff(G.cells.facePos)==4);
plotGrid(G,'FaceColor','none','EdgeAlpha',.05);
plotGrid(G,id,'FaceColor','r','EdgeAlpha',.5,'MarkerSize',5,'Marker','o');




%% saigup model

grdecl = fullfile(getDatasetPath('SAIGUP'), 'SAIGUP.GRDECL');
grdecl = readGRDECL(grdecl);

% Create the grid data structure.
G = processGRDECL(grdecl,'PreserveCpNodes',true,'Tolerance',0);

% Plot the resulting geometry.
% figure (1), hg1 = plotGrid(G, 'FaceAlpha', 0.625);
% view(3), grid on, axis tight

plotGrid(G,'FaceColor','none','EdgeAlpha',.1);

plotFaces(G, G.faces.tag>0, ...
'FaceColor','red','FaceAlpha',.2, ...
'EdgeColor','r','EdgeAlpha',.1);
axis off; view(-155,80); zoom(1.7);

%%

grdecl = pinchedLayersGrdecl([50 50 10], @(x) .05+0.07*x);

% Create the grid data structure
G = computeGeometry(processGRDECL(grdecl));

% Plot the geometry
hg = plotGrid(G, 'FaceAlpha', 0.625);
view(3), grid on, axis tight


%% 

% 1) Create a 20-by-20-by-5 faulted grid with a fault drop of 0.15 (m).
grdecl = simpleGrdecl([6, 5, 100], 0.15);

msh = GRDECLtoGRES(grdecl,0.01);

% Create the grid data structure.
G = processGRDECL(grdecl,'PreserveCpNodes',true,'Tolerance',0.01);

% Plot the resulting geometry.
figure (1), hg1 = plotGrid(G, 'FaceAlpha', 0.625);
view(3), grid on, axis tight