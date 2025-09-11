% Tolerance error between the results.
errorTol = 1e-2;

% Get the full path of this test
testPath = mfilename('fullpath');
cd(fileparts(testPath));

model = ModelType("VariabSatFlow_FVTPFA");
simParam = SimulationParameters(fullfile('simParam.xml'),model);

% Create the topology object
topology = Mesh();
topology.importGMSHmesh(fullfile('Mesh','Column.msh'));
elems = Elements(topology,2);
faces = Faces(model,topology);
grid = struct('topology',topology,'cells',elems,'faces',faces);

% Degree of freedom manager
dofmanager = DoFManager(topology,model);

% Creating boundaries conditions.
bound = Boundaries("BC.dat",model,grid);

% to set initial condition.
z = elems.mesh.cellCentroid(:,3);
wLev = 9.;

sol = struct('time', [], 'pressure', [],'saturation', []);
listMat = ["matTab.dat" "matAna.dat"];
for sim = 1:numel(listMat)
  run("runRichards.m");
  sol(sim).time = [printUtils.results.expTime]';
  sol(sim).time = sol(sim).time(2:end);
  sol(sim).pressure = [printUtils.results.expPress]';
  sol(sim).saturation = [printUtils.results.expSat]';

  clearvars mat domain Solver simState printUtils
end

clearvars -except sol errorTol
for sim = 1:numel(sol(1).time)
  % Check if the solution are evaluated at the same time.
  err = norm((sol(1).time(sim)-sol(2).time(sim))./sol(2).time(sim));
  assert( err < errorTol, ['Test - %i - Results are not evaluated at' ...
    ' the same time - tabular(%i) analitical(%i)\n'], sim, ...
    sol(1).time(sim), sol(2).time(sim));

  % Check if the pressure are approximately the same.
  err = norm((sol(1).pressure(sim,:)-sol(2).pressure(sim,:))./sol(2).pressure(sim,:));
  assert( err < errorTol, ['Test - %i - Results of pressure are not' ...
    ' the same - tabular(|%i|) analitical(|%i|)\n'], sim, ...
    norm(sol(1).pressure(sim,:)), norm(sol(2).pressure(sim,:)));

  % Check if the saturation are approximately the same.
  err = norm((sol(1).saturation(sim,:)-sol(2).saturation(sim,:))./sol(2).saturation(sim,:));
  assert( err < errorTol, ['Test- %i - Results of saturation are not the' ...
    ' same - tabular(|%i|) analitical(|%i|)\n'], sim, ...
    norm(sol(1).saturation(sim,:)),norm(sol(2).saturation(sim,:)));
end