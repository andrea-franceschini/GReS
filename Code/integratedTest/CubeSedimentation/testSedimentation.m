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
ref(1).press   = [5.696736937229210e-04; 0];
ref(1).stress  = [-1.650006000161258e+03; -1.499882008974347e+02];
ref(1).strain  = [-3.762677543472629e-06; 0];

ref(2).time = 10.1;
ref(2).press   = [1.392902213852669e-03; 7.974303343768915e-04];
ref(2).stress  = [-3.014897805099393e+03; -1.514880031633756e+03];
ref(2).strain  = [-2.736522102337879e-05; 0];

ref(3).time = 11.;
ref(3).press   = [1.464902159932661e-03; 8.680868245770957e-04; 0];
ref(3).stress  = [-3.149887113907138e+03; -1.649869341784957e+03; -1.499691689054160e+02];
ref(3).strain  = [-2.929151164702439e-05; 0; 0];

ref(4).time = 21.;
ref(4).press   = [2.475626627689976e-03; 2.071332546363613e-03; 1.070227202550541e-03; 0];
ref(4).stress  = [-4.649768112157017e+03; -3.149750147513581e+03; -1.649850107652560e+03; -1.499572234211735e+02];
ref(4).strain  = [-4.466828432299391e-05; -2.552886492208846e-05; 0; 0];

ref(5).time = 30.;
ref(5).press   = [3.452834660752819e-03; 3.139508602098363e-03; 2.408416532527968e-03; 1.118558686820950e-03];
ref(5).stress  = [-5.999660943025895e+03; -4.499642887414438e+03; -2.999742577540141e+03; -1.499849912939399e+03];
ref(5).strain  = [-5.373249816442258e-05; -3.827222210211486e-05; -2.157771238235775e-05; 0];

% Get the full path of this test
input_dir = 'Input/';
testPath = mfilename('fullpath');
cd(fileparts(testPath));

simParam = SimulationParameters('Start',0.,'End',30.1,'DtInit',1e-0,...
  'DtMin',1e-0,'DtMax',1e2,'incrementFactor',1.0);

% Create an object of the Materials class and read the materials file
mat = Materials(fullfile(input_dir,'materials.xml'));

% Create and set the print utility
printUtils = OutState(fullfile(input_dir,'output.xml'));

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

