% Result's used as reference to comparison.
ref = repelem(struct('press', 1, 'stress', 1, 'strain', 1), 3);

ref(1).press   = 3.975934704651128e-04;
ref(1).stress  = -2.514982394451165e+03;
ref(1).strain  = -3.640572790216300e-05;

ref(2).press   =      [0.001154246105681; 0.000667776770957];
ref(2).stress  = 1e+3*[-4.014964600218962; -2.364983828025822];
ref(2).strain  = 1e-4*[-0.548771014813824; -0.339770784178433];

ref(3).press   =      [0.001986976678321; 0.001645170303570; 0.000815055888433];
ref(3).stress  = 1e+3*[-5.514946729908823; -3.864965813052720; -2.364983680746704];
ref(3).strain  = 1e-4*[-0.674120817397354; -0.533762716145965; -0.339839911098807];

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
solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils);
solver.simulationLoop();

% Solution
sol = printUtils.results;
time = printUtils.timeList;
clearvars -except sol errorTol time ref

% Tolerance error between the results.
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

  % Check the strain result.
  anaVal = sol(sim).strain(1:4:end);
  anaRef = ref(sim).strain;
  err = norm((anaVal-anaRef)./anaRef);
  assert( err < errorTol, ['Simulation at Time - %i \n' ...
    'Strain Result   - %s\n' ...
    'Strain Expected - %s\n' ...
    'Error - %d\n'], sol(sim).time, mat2str(anaVal), mat2str(anaRef), err);
end