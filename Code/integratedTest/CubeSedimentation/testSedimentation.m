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
ref(1).press   = [5.662777838171824e-05; 0];
ref(1).stress  = [1.659881522908224e+02; 1.015000008021166e+03];
ref(1).strain  = [-3.739167694651991e-06; 0];

ref(2).time = 10.1;
ref(2).press   = [9.438560935391386e-05; 3.503567426417822e-05];
ref(2).stress  = [3.024881875256033e+02; 1.151500045978104e+03];
ref(2).strain  = [-2.723513899440930e-05; 0];

ref(3).time = 11.;
ref(3).press   = [9.730716827343308e-05; 3.781508100371161e-05; 0];
ref(3).stress  = [3.159881918230939e+02; 1.165000050417747e+03; 1.014998112222454e+03];
ref(3).strain  = [-2.915516776363712e-05; 0; 0];

ref(4).time = 21.;
ref(4).press   = [1.303685685634306e-04; 9.002697772995743e-05; 3.538796715886759e-05; 0];
ref(4).stress  = [4.659882389733550e+02; 1.315000078417512e+03; 1.164998157046149e+03; 1.014996920553215e+03];
ref(4).strain  = [-4.449125071534490e-05; -4.781830003829677e-06; 0; 0];

ref(5).time = 30.;
ref(5).press   = [1.643366570224716e-04; 1.330572243731780e-04; 8.881109714649322e-05; 3.303886677707911e-05];
ref(5).stress  = [6.009882771957615e+02; 1.450000107577761e+03; 1.299998175813514e+03; 1.149996959704843e+03];
ref(5).strain  = [-5.353770183370028e-05; -8.229695816906608e-06; -3.870681063557507e-06; 0];

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

