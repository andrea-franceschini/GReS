% model = ModelManager();
% model.createModel('terzaghi.xml');
% 
% 
% model.runProblem();


% new model without requiring ModelType and DoFManager. Input file for physicsSolver
% is required

simparam = SimulationParameters('simparam.xml');

mesh = Mesh();
mesh.importMesh('Mesh/Column_hexa.msh');

elems = Elements(mesh,2);
faces = Faces(mesh);

grid = struct('topology',mesh,'cells',elems,'faces',faces);

mat = Materials('materials.xml');

bc = Boundaries('boundaryConditions.xml',grid);

printUtils = OutState(mesh,'output.xml');


% create the Discretizer (key-value pair)
domain = DiscretizerNew('grid',grid,...
                        'materials',mat,...
                        'boundaries',bc,...
                        'outstate',printUtils);

domain.addPhysicsSolver('solver.xml');

solv = FCSolver(simparams,domain);

stat = solv.NonLinearLoop();

domain.outstate.finalize();
