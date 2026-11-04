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
linearSolver.solveSingle('matrices/stickslipEFEM.mat');



% twoPhase


% multiPhasePoromechanics





