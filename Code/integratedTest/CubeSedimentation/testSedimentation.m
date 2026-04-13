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
ref(1).press   = [ 1.120006181306114e-03; 0];
ref(1).stress  = [-2.499981842414250e+03;-1.000000000000000e+03];
ref(1).strain  = [-3.617566676766528e-05; 0];

ref(2).time = 10.1;
ref(2).press   = [ 1.196159396145166e-03; 1.172110608733250e-04];
ref(2).stress  = [-2.514981595885240e+03;-1.014999712413143e+03];
ref(2).strain  = [-3.640571542542461e-05;-5.517824079168186e-07];

ref(3).time = 20.;
ref(3).press   = [ 2.291947815971757e-03; 1.803647813734837e-03; 0];
ref(3).stress  = [-3.999963632893046e+03;-2.499981158772618e+03;-1.000000000000000e+03];
ref(3).strain  = [-5.473174585603535e-05;-3.617565597127934e-05; 0];

ref(4).time = 30.;
ref(4).press   = [ 3.486696597097168e-03; 3.135308537769740e-03; 2.295614702488390e-03];
ref(4).stress  = [-5.499945400564700e+03;-3.999962789532326e+03;-2.499980666805729e+03];
ref(4).strain  = [-6.730454262555916e-05;-5.473173753180543e-05;-3.617564820190809e-05];

ref(5).time = 30.;
ref(5).press   = [ 3.486696597097168e-03; 3.135308537769740e-03; 2.295614702488390e-03; 0];
ref(5).stress  = [-5.499945400564700e+03;-3.999962789532326e+03;-2.499980666805729e+03;-1.000000000000000e+03];
ref(5).strain  = [-6.730454262555916e-05;-5.473173753180543e-05;-3.617564820190809e-05; 0];

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

