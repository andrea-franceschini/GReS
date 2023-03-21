clc;
close all;
clear;

addpath('../read');
mesh = Mesh();
mesh.importGMSHmesh('../read/mesh.msh');
%mxVTKWriter(mesh.m_coordinates, mesh.m_cells, mesh.m_cellVTKType, mesh.m_cellNumVerts);
mesh.finalize();
pointData3D.name = 'val1';
pointData3D.data = mesh.coordinates(:,2);
pointData2D.name = 'val2';
pointData2D.data = mesh.coordinates(:,3);
cellData3D = repmat(struct('name', 1, 'data', 1), 2, 1);
cellData3D(1).name = 'val3';
cellData3D(1).data = mesh.cellCentroid(:,3);
cellData3D(2).name = 'mat';
cellData3D(2).data = [1 : length(mesh.cellCentroid(:,3))]';
cellData2D.name = 'val4';
cellData2D.data = mesh.surfaceCentroid(mesh.findSurfacesOfRegion('f1'),:);

V = VTKOutput(mesh);
V.setSurfaces({'f1'});
V.writeVTKFile(2.3, pointData3D, cellData3D, pointData2D, cellData2D);
