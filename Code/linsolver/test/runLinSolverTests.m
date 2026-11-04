% spflow
solveSingle('matrices/spflow_fluxBarrier.mat');

%% mechanics 
solveSingle('matrices/mech_treporti.mat');


%% richards
solveSingle('matrices/richards_column.mat');

%% biot (surfaceLoad)
solveSingle('matrices/biotFEM_treporti.mat');


%% stick slip (simulation, EW, linear)
run('simul/runStickSlipOpen.m')

% EFEM
solveSingle('matrices/stickslipEFEM.mat');



% twoPhase


% multiPhasePoromechanics





