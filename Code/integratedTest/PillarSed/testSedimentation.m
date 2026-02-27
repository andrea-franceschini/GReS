input_dir = 'Input/';

% Get the full path of this test
testPath = mfilename('fullpath');
cd(fileparts(testPath));

simParam = SimulationParameters(fullfile(input_dir,'simparam.xml'));

% Create an object of the Materials class and read the materials file
mat = Materials(fullfile(input_dir,'materials.xml'));

% Create and set the print utility
printUtils = OutState(fullfile(input_dir,'output.xml'));

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Materials',mat);
domain.addPhysicsSolver(fullfile(input_dir,'solver.xml'));

% Solve the problem
solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils);
solver.simulationLoop();

% Solution
sol = printUtils.matFile;
time = printUtils.timeList;
clearvars -except sol errorTol time

% Tolerance error between the results.
ref = repelem(struct('press', 1, 'stress', 1, 'compact', 1), 5);

ref(1).press   = 0.001942180634281;
ref(1).stress  = -3.024998057819366e3;
ref(1).compact = -0.115585318301551e-5;

ref(2).press   =      [0.004118153264475; 0.002495189371988];
ref(2).stress  = 1e+3*[-7.389995881846737; -3.639997504810628];
ref(2).compact = 1e-4*[-0.214558042919429; -0.268161428382395];

ref(3).press   =      [0.003101270267381; 0.001955474882098; 0.000336751341030];
ref(3).stress  = 1e+4*[-1.038999689872973; -0.663999804452512; -0.288999966324866];
ref(3).compact = 1e-4*[-0.291995625179768; -0.482231001934589; -0.483386844408959];

ref(4).press   =      [0.009580919101506; 0.008783840785885; 0.006917637394896; 0.003428371777083];
ref(4).stress  = 1e+4*[-1.488999041908090; -1.113999121615922; -0.738999308236260; -0.363999657162822];
ref(4).compact = 1e-4*[-0.373779802171336; -0.681616130293312; -0.896174087113224; -0.949777414299979];

ref(5).press   = 1e-5*[0.400107850725799; 0.368673964675891; 0.301329837932782; 0.196735583165046; 0.065574214007363];
ref(5).stress  = 1e+4*[-1.787499999599892; -1.412499999631326; -1.037499999698670; -0.662499999803264; -0.287499999934426];
ref(5).compact = 1e-3*[-0.041530671751304; -0.077710053007387; -0.106876998892472; -0.125849659246620; -0.125849659241437];

errorTol = 1e-8;
for sim = 1:numel(sol)
  numcell = numel(sol(sim).pressure);

  % Check if the growing is uniform
  assert(mod(numcell,4)==0,'The sedimentation was not uniform in the test.');

  % Check the pressure result.
  anaVal = sol(sim).pressure(1:4:end);
  anaRef = ref(sim).press;
  err = norm((anaVal-anaRef)./anaRef);
  assert( err < errorTol, ['Simulation at Time - %i \n' ...
    'Pressure Result   - %s\n' ...
    'Pressure Expected - %s\n' ...
    'Error - %d\n'], sol(sim).time, mat2str(anaVal), mat2str(anaRef), err);  

  % Check the stress result.
  anaVal = sol(sim).stress(1:4:end);
  anaRef = ref(sim).stress;
  err = norm((anaVal-anaRef)./anaRef);
  assert( err < errorTol, ['Simulation at Time - %i \n' ...
    'Stress Result   - %s\n' ...
    'Stress Expected - %s\n' ...
    'Error - %d\n'], sol(sim).time, mat2str(anaVal), mat2str(anaRef), err);

  % Check the compaction result.
  anaVal = sol(sim).compaction(1:4:end);
  anaRef = ref(sim).compact;
  err = norm((anaVal-anaRef)./anaRef);
  assert( err < errorTol, ['Simulation at Time - %i \n' ...
    'Compaction Result   - %s\n' ...
    'Compaction Expected - %s\n' ...
    'Error - %d\n'], sol(sim).time, mat2str(anaVal), mat2str(anaRef), err);
end