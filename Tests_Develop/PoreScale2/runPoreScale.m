clear
close all
clc

%%
runFluidProblem;
press = solverFlow.getState().data.pressure;
clearvars -except press solverFlow gridFluid

%%
fluid2grains;

clearvars -except press_surf press_grain mshGrain press_grain_inner


plotFunction(mshGrain,'OUT_press_grain',press_grain_inner)



%%
runGrainsProblem;