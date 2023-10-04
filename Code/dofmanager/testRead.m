clear; clc; close all;
topology = Mesh();
%setting the model and physics included
model = ModelType(["SinglePhaseFlow_FEM","Poromechanics_FEM"]);
% Set the input file name

%%%% SET UP THE TEST GRID

fileName = 'TestDoFManager.msh';
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);
elems = Elements(topology);
faces = Faces(model, topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);

%%%TESTING DOF MANAGER CLASS
tic
fileName = 'dof.dat';
dofmanager = DoFManager(topology, model, fileName);
dof_time = toc;
test = dofmanager.getDofTables();
testdof = dofmanager.getLocDoF('Flow');
%tab = getSubTable(dofmanager,2);

%%% TEST BCs ACCESS TO DOFS 
% Set the input file
fileName = ["dir_BC_flow_tetra.dat","dir_BC_poro_tetra.dat","neuSurf_BC_poro_tetra.dat"];
%
% Create an object of the "Boundaries" class and read the boundary
% conditions
bound = Boundaries(fileName,model,grid,dofmanager);

keys = bound.db.keys;
for i=1:length(keys)
    testDOF{i} = bound.getDofs(keys{i});
    testENTS{i} = bound.getEntities(keys{i});
end

