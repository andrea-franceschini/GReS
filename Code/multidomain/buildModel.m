function varargout = buildModel(fileName)

% build a GReS simulation from a unique input file


outStruct = readstruct(fileName,AttributeSuffix="");

if isfield(outStruct,"Domain")
  outStruct = outStruct.Domain;
end

nD = numel(outStruct);
domains = [];

for i = 1:nD
  domains = [domains; defineDomain(outStruct(i))];
end

varargout{1} = domains;

if nargout > 1
interfaces = buildInterfaces(fileName,domains);
varargout{2} = interfaces;
end

end


function domain = defineDomain(input)

topology = Mesh();
geom = input.Geometry;
meshFile = getXMLData(geom,[],"fileName");
topology.importMesh(meshFile);

if isfield(input,"Materials")
  mat = Materials(input);
else
  mat = [];
end

gNPoints = 0;
if isfield(input,"Gauss")
  g = input.Gauss;
  gNPoints = getXMLData(g,[],'nGP');
end

if gNPoints~=0
  elems = Elements(topology,gNPoints);
else
  elems = Elements(topology);
end

faces = Faces(topology);
grid = struct('topology',topology,'cells',elems,'faces',faces);


% boundary conditions
if isfield(input,'BoundaryConditions')
  bound = Boundaries(input,grid);
else
  bound = [];
end

% output manager
if isfield(input,"Output")
  printUtils = OutState(grid.topology,input);
else
  printUtils = OutState(grid.topology);
end


domain = Discretizer('grid',grid,...
                     'materials',mat,...
                     'boundaries',bound,...
                     'outstate',printUtils);

domain.addPhysicsSolver(input.Solver);


end

