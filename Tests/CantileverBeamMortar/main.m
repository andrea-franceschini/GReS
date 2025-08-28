mshMaster = Mesh();
mshSlave = Mesh();
mshMaster.importMesh('Mesh/leftBeam.vtk');
mshSlave.importMesh('Mesh/right_beam.vtk');

interfMaster = getSurfaceMesh(mshMaster,2);
interfSlave = getSurfaceMesh(mshSlave,1);

elM = Elements(interfMaster,2);
elS = Elements(interfSlave,2);

m = Mortar(interfMaster,interfSlave);