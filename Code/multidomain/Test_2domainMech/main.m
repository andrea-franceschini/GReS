clear
close all

%% Testing multidomain model and inherent novel methods
fileName = "simParam.dat";
model = ModelType("Poromechanics_FEM");
simParam = SimulationParameters(model,fileName);

fileName = "domains.dat";
% All the preprocessing is performed on different domain
% All istances of the classes are collected in a structure array
modelStructure = buildModelStruct(model,simParam,fileName);

% SOLUTION ALGORITHM



