% spflow
linearSolver.solveSingle('matrices/spflow_fluxBarrier.mat');

%% mechanics 
linearSolver.solveSingle('matrices/mech_treporti.mat');


%% richards
linearSolver.solveSingle('matrices/richards_column.mat');

%% biot (surfaceLoad)
linearSolver.solveSingle('matrices/biotFEM_treporti.mat');


%% stick slip (simulation, EW, linear)
run('simul/runStickSlipOpen.m')
linearSolver.solveSingle('matrices/stickslipmat.mat');

%% EFEM
[~,out] = linearSolver.solveSingle('matrices/stickslipEFEM.mat');

%% twoPhase
linearSolver.solveSingle('matrices/twoPhase_SPE.mat');

%% multiPhasePoromechanics
linearSolver.solveSingle('matrices/mpPoro_staircase.mat');




