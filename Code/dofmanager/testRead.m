clear; clc; close all;
topology = Mesh();
%
% Set the input file name
fileName = 'TestDoFManager.msh';
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);

fileName = 'dof.dat';
dofmanager = DoFManager(topology,fileName);

tab = getSubTable(dofmanager,2);