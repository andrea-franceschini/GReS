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
rmpath(genpath('Output'));

% Result's used as reference to comparison.
ref = repelem(struct('press', 1, 'stress', 1, 'strain', 1), 3);

ref(1).time = 10.;
ref(1).press   = 1e-3*[ 0.399993514141854; 0.];
ref(1).stress  = 1e+3*[-2.499982562426917; -1.000000000000000];
ref(1).strain  = 1e-4*[-0.361756781384363; 0.];

ref(2).time = 10.1;
ref(2).press   = 1e-3*[ 0.548147483547348; 0.117210890155328];
ref(2).stress  = 1e+3*[-2.514982243897152; -1.014999712413314];
ref(2).strain  = 1e-4*[-0.364057256591261; -0.005517824079227];

ref(3).time = 20.;
ref(3).press   =      [ 0.001132446099120; 0.000644146008496; 0.];
ref(3).stress  = 1e+3*[-3.999964792394763;-2.499982318274423;-1.000000000000000];
ref(3).strain  = 1e-4*[-0.547317573006713;-0.361756742826682; 0.];

ref(4).time = 30.;
ref(4).press   =      [ 0.002010927083071; 0.001659538887333; 0.000819844775605];
ref(4).stress  = 1e+3*[-5.499946876334215;-3.999964265301975;-2.499982142575656];
ref(4).strain  = 1e-4*[-0.673045532192377;-0.547317520981050;-0.361756715079526];

ref(5).time = 30.;
ref(5).press   =      [ 0.002010927083071; 0.001659538887333; 0.000819844775605; 0];
ref(5).stress  = 1e+3*[-5.499946876334215;-3.999964265301975;-2.499982142575656;-1.000000000000000];
ref(5).strain  = 1e-4*[-0.673045532192377;-0.547317520981050;-0.361756715079526; 0];

% Get the full path of this test
input_dir = 'Input/';
testPath = mfilename('fullpath');
cd(fileparts(testPath));

simParam = SimulationParameters(fullfile(input_dir,'simparam.xml'));

% Create an object of the Materials class and read the materials file
mat = Materials(fullfile(input_dir,'materials.xml'));

% Create and set the print utility
printUtils = OutState(fullfile(input_dir,'output.xml'));

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Materials',mat);
domain.addPhysicsSolvers(fullfile(input_dir,'solver.xml'));

% Solve the problem
solver = EvolvingGrid('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils, ...
                           'growprint',true);
solver.simulationLoop();

% Solution
sol = printUtils.results;
time = printUtils.timeList;
clearvars -except sol errorTol time ref

% Tolerance error between the results.
errorTol = 1e-8;
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

