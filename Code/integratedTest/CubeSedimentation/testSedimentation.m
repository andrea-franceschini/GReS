% This unit test evaluates the sedimentation of a region composed of four
% cells. The formation of two new layers due to sediment deposition
% is expected.
% The simulation is configured to output results whenever the mesh grows,
% as well as at two specific times (one immediately after a new layer
% forms and another exactly at the moment a new layer is generated).
% Therefore, an output with five results is expected. However, for the
% specific time corresponding exactly to the formation of a new layer,
% due to the internal logic of GRES, this result will be recorded before
% the creation of the new layer.

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
scriptDir = fileparts(scriptFullPath);
cd(scriptDir);
rmpath(genpath(fullfile(gres_root,...
        "/Code/integratedTest/CubeSedimentation/Output")));
% Result's used as reference to comparison.
ref = repelem(struct('press', 1, 'stress', 1, 'strain', 1), 3);

ref(1).time = 1.;
ref(1).press   = [5.6971e-05; 0];
ref(1).stress  = [1.6499e+02; 1.5000e+01];
ref(1).strain  = [-3.7629e-06; 0];

ref(2).time = 10.1;
ref(2).press   = [1.3929e-04; 7.9747e-05];
ref(2).stress  = [3.0149e+02; 1.5150e+02];
ref(2).strain  = [-2.7366e-05; 0];

ref(3).time = 11.;
ref(3).press   = [1.4649e-04; 8.6813e-05; 0];
ref(3).stress  = [3.1499e+02; 1.6500e+02; 1.4998e+01];
ref(3).strain  = [-2.9292e-05; 0; 0];

ref(4).time = 21.;
ref(4).press   = [2.4757e-04; 2.0714e-04; 1.0703e-04; 0];
ref(4).stress  = [4.6499e+02; 3.1500e+02; 1.6500e+02; 1.4997e+01];
ref(4).strain  = [-4.4669e-05; -2.5531e-05; 0; 0];

ref(5).time = 30.;
ref(5).press   = [3.4530e-04; 3.1397e-04; 2.4086e-04; 1.1186e-04];
ref(5).stress  = [5.9999e+02; 4.5000e+02; 3.0000e+02; 1.5000e+02];
ref(5).strain  = [-5.3733e-05; -3.8275e-05; -2.1579e-05; 0];

% Get the full path of this test
input_dir = 'Input/';
testPath = mfilename('fullpath');
cd(fileparts(testPath));

simParam = SimulationParameters('Start',0.,'End',30.1,'DtInit',1e-0,...
  'DtMin',1e-0,'DtMax',1e2,'incrementFactor',1.0);

% Create an object of the Materials class and read the materials file
mat = Materials(fullfile(input_dir,'materials.xml'));

% Create and set the print utility
% printUtils = OutState(fullfile(input_dir,'output.xml'));
printUtils = OutState('outputFile','Output/Results','printTimes',[10.1,30]);

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('materials',mat);
domain.addPhysicsSolvers(fullfile(input_dir,'solver.xml'));

% Solve the problem
solver = EvolvingGrid('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils, ...
                           'growprint',1);
solver.simulationLoop();

% Solution
sol = printUtils.results;
time = printUtils.timeList;
clearvars -except sol errorTol time ref

% Tolerance error between the results.
errorTol = 1e-4;
epsVal = 1e-12;
for sim = 1:numel(sol)
  numcell = numel(sol(sim).pressure);

  % Check if the growing is uniform
  assert(mod(numcell,4)==0,'The sedimentation was not uniform in the test.');

  % Check the same time.
  anaVal = sol(sim).time;
  anaRef = ref(sim).time;
  assert( anaVal == anaRef, ['Simulation not at the same Time - ' ...
    'Simulation Time - %s\n' ...
    'Expected   Time - %s\n'], anaVal, anaRef);  

  % Check the pressure result.
  anaVal = sol(sim).pressure(1:4:end);
  anaRef = ref(sim).press;  
  err = norm((anaVal-anaRef)./(anaRef + epsVal));
  assert( err < errorTol, ['Simulation at Time - %i \n' ...
    'Pressure Result   - %s\n' ...
    'Pressure Expected - %s\n' ...
    'Error - %d\n'], sol(sim).time, mat2str(anaVal), mat2str(anaRef), err);  

  % Check the stress result.
  anaVal = sol(sim).stress(1:4:end);
  anaRef = ref(sim).stress;
  err = norm((anaVal-anaRef)./(anaRef + epsVal));
  assert( err < errorTol, ['Simulation at Time - %i \n' ...
    'Stress Result   - %s\n' ...
    'Stress Expected - %s\n' ...
    'Error - %d\n'], sol(sim).time, mat2str(anaVal), mat2str(anaRef), err);

  % Check the strain result.
  anaVal = sol(sim).strain(1:4:end);
  anaRef = ref(sim).strain;
  err = norm((anaVal-anaRef)./(anaRef + epsVal));
  assert( err < errorTol, ['Simulation at Time - %i \n' ...
    'Strain Result   - %s\n' ...
    'Strain Expected - %s\n' ...
    'Error - %d\n'], sol(sim).time, mat2str(anaVal), mat2str(anaRef), err);
end

s = readstruct("Output/Results.pvd","FileType","xml");
v = s.Collection.DataSet;
assert( length(v) == 5, 'Five results is expected for this Test!');