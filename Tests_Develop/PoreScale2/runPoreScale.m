clear
close all
clc

%%
runFluidProblem;
press = solverFlow.getState().data.pressure;
clearvars -except press solverFlow gridFluid

%%
fluid2grains;


plotFunction(mshGrain,'OUT_press_grain',press_grain)



%%
runGrainsProblem;