%
clear
close all

grid = "simple"; % simple, norne 
outName = "gridTest";
%

NX = 60; 
NY = 40;
NZ = 10;

switch grid
  case 'simple'

    %grdecl = pinchedLayersGrdecl([20, 15, 10], 0.08);
    %grdecl = threeLayers(10, 10, [5, 10, 5]);
    %grdecl = makeModel3([100, 60, 15]);
    %grdecl = oneSlopingFault([90, 10, 16], 5);
    grdecl = simpleGrdecl([NX,NY,NZ], 0.025);
  case 'norne'
    grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');
    grdecl = readGRDECL(grdecl);
    usys   = getUnitSystem('METRIC');
    grdecl = convertInputUnits(grdecl, usys);

    % Save ACTNUM until later and then override the ACTNUM values to obtain a
    % grid that contains all cells
    actnum        = grdecl.ACTNUM;
    grdecl.ACTNUM = ones(prod(grdecl.cartDims),1);
end

%% testing corner point grid modification
grdecl = refineLayers(grdecl,5,5);
dir = 1;
target = round(NX/2)+1;
geom_factor = 0.93;
grdecl = movePillars(grdecl, dir , target, 1.05, "before");
grdecl = movePillars(grdecl, dir , target, geom_factor, "after");
grdecl = movePillars2(grdecl, 2 , [0 0.5 1], [0.2, 0.2 0.2], 0.005, 0.04);
grdecl = removePillars(grdecl, dir, 1:(NX/2-4));
%%


G = processGRDECL(grdecl,'PreserveCpNodes',true,'Tolerance',0);
G = computeGeometry(G,'CpGeometry',true);


%%

[mesh,conn] = GRDECLtoGRES(grdecl,G,0);

% rescale the base grid to 1x1x0.5
scale = [1e3,1e3,100];
mesh.coordinates = mesh.coordinates.*scale;

%%
elem = Elements(mesh,4);

%%
plotFunction(mesh,outName,ones(mesh.nNodes,1));



mesh.nSurfaceTag = max(mesh.surfaceTag);
%%
mshSlave = getSurfaceMesh(mesh,2);
mshMaster = getSurfaceMesh(mesh,1);

plotFunction(mshSlave,outName+"slave",ones(mshSlave.nNodes,1));

plotFunction(mshMaster,outName + "master",ones(mshMaster.nNodes,1));
