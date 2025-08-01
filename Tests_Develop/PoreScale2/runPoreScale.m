clear
close all
clc

%%
runFluidProblem;
press = solverFlow.getState().data.pressure;
clearvars -except press gridFluid linSystFluid dofFluid simParam

%%
fluid2grains;

runGrainsProblem;