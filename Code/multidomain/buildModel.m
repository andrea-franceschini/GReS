function varargout = buildModel(fileName)
% This function constructs a full simulation model from a single XML input
% file. The goal is to assemble all components required for a GReS
% simulation (geometry, mesh, materials, integration rules, boundary
% conditions, and output settings) into an array of Discretizer object, one
% for each domain in the simulation.
% Optional interface construction is also supported.
% 
% FUNCTION: buildModel(fileName)
% 
% Reads the XML input using readstruct, extracts the <Domain> block when
% present, and iterates over each domain specification. For each domain
% entry, defineDomain is called to create a Discretizer object.
% The function returns:
% - The list of fully constructed domain objects  
% - Optionally, the interface objects connecting multiple domains, if any.
% 
% Example usage:
% domains = buildModel("input.xml");
% [domains, interfaces] = buildModel("input.xml");



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

if ~any(topology.cellVTKType==10)
  faces = Faces(topology);
else
  faces = [];
end

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

