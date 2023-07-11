close all
clear
addpath('../modeltype');
% model = ModelType("VariabSatFlow_FVTPFA");
model = ModelType("SinglePhaseFlow_FVTPFA");
fileName = "SimParamSPF.dat";
simParam = SimulationParameters(model,fileName)