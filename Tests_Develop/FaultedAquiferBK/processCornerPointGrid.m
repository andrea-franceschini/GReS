function [mesh,newDims] = processCornerPointGrid(outName,dims,numbRockCells,scale)


NX = dims(1); NY = dims(2); NZ = dims(3);
% layer_ref = ref(1);
% n_refs = ref(2);

grdecl = simpleGrdecl([NX,NY,NZ], 0.01);

% modifying original MRST corner point grid
%grdecl = refineLayers(grdecl,layer_ref,n_refs);
dir = 1;
target = round(NX/2)+1;
geom_factor = 0.85;
grdecl = movePillars(grdecl, dir , target, 1.03, "before");
grdecl = movePillars(grdecl, dir , target, geom_factor, "after");
% grdecl = movePillars2(grdecl, 2 , [0 0.5 1], [0.2, 0.2 0.2], 0.005, 0.04);
grdecl = movePillars2(grdecl, 2 , [0 1], [0.5, 0.5], 0.01, 0.1);
grdecl = removePillars(grdecl, dir, 1:(NX/2-numbRockCells));

G = processGRDECL(grdecl,'PreserveCpNodes',true,'Tolerance',0);
G = computeGeometry(G,'CpGeometry',true);

newDims = G.cartDims;

[mesh,~] = GRDECLtoGRES(grdecl,G,0);

% rescale the base grid to 1x1x0.5
mesh.coordinates = mesh.coordinates.*scale;

plotFunction(mesh,outName,ones(mesh.nNodes,1));

mesh.nSurfaceTag = max(mesh.surfaceTag);
%%
mshSlave = getSurfaceMesh(mesh,2);
mshMaster = getSurfaceMesh(mesh,1);

if ~isempty(mshSlave)

% spot if there are inactive surfaces 
if mshSlave.nSurfaces < newDims(2)*newDims(3)
  error('Mortar inactive surfaces in the fault. Consider reducing the fault jump')
end

plotFunction(mshSlave,outName+"slave",ones(mshSlave.nNodes,1));
plotFunction(mshMaster,outName + "master",ones(mshMaster.nNodes,1));

end

end

