clear; clc; close all;
topology = Mesh();
%setting the model and physics included
model = ModelType(["SinglePhaseFlow_FVTPFA","Poromechanics_FEM"]);
% Set the input file name
fileName = 'TestDoFManager.msh';
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);

fileName = 'dof.dat';
dofmanager = DoFManager(topology,fileName,model);
test = dofmanager.getDofTables();
testdof = dofmanager.getDoF([1 2 6 9],'Flow');
%tab = getSubTable(dofmanager,2);