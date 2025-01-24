close all;
clear;

% -------------------------- SET THE PHYSICS -------------------------
model = ModelType(["SinglePhaseFlow_FVTPFA","Poromechanics_FEM"]);
%
% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = "simParam.dat";
simParam = SimulationParameters(fileName);
%
% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();
%
% Set the input file name
fileName = 'ReservoirTest_Hexa.msh';
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);

%
% Degree of freedom manager 
fname = 'dof.dat';
dofmanager = DoFManager_new(topology,model,fname);

% Testing DoF manager
out1 = dofmanager.getDoF("Poromechanics");
out2 = dofmanager.getDoF("SPFlow");
cTags = dofmanager.getFieldTags("Poromechanics");