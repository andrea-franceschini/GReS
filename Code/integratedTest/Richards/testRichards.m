% Tolerance error between the results.
errorTol = 1e-2;

% Get the full path of this test
testPath = mfilename('fullpath');
cd(fileparts(testPath));


listMat = ["tabular","analytical"];

s = struct('time', [], 'pressure', [],'saturation', []);

for i = 1:numel(listMat)
  s(i) = run(listMat(i)); 
end

validate(s);



function sol = run(curveType)
% model = ModelType("VariabSatFlow_FVTPFA");
simParam = SimulationParameters(fullfile('Input','simParam.xml'));

% Create the topology object
topology = Mesh();
topology.importMesh(fullfile('Input','Mesh','Column.msh'));
elems = Elements(topology,2);
faces = Faces(topology);
grid = struct('topology',topology,'cells',elems,'faces',faces);

% Creating boundaries conditions.
bound = Boundaries(grid,"Input/richardsBCs.xml");

% to set initial condition.
z = elems.mesh.cellCentroid(:,3);
wLev = 9.;

mat = Materials('Input/mat.xml');

matName = mat.getMaterialNames;
switch curveType
  case 'tabular'
    relPath = "Code/integratedTest/Richards/Input/Materials/krCurve_2000.dat";
    capPath = "Code/integratedTest/Richards/Input/Materials/pcCurve_2000.dat";
    mat.addCapillaryCurves(matName,'type',"tabular",...
      'relativePermeabilityPath',relPath,...
      'capillaryCurvePath',capPath);
  case 'analytical'
    mat.addCapillaryCurves(matName,'type',"preset",'soilName',"sand");
end


% Create and set the print utility
printUtils = OutState('printTimes',[5;10]);


% Create object handling construction of Jacobian and rhs of the model
% linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);
domain = Discretizer('Grid',grid,...
                     'Materials',mat,...
                     'Boundaries',bound);

domain.addPhysicsSolver('VariablySaturatedFlow');

% set initial conditions directly modifying the state object
domain.state.data.pressure = getFluid(mat).getSpecificWeight()*(wLev-z);

% Solve the problem
solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils);
solver.simulationLoop();


sol = solver.output.results;

end




function validate(sol)

assert(numel(sol)==2);

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

end