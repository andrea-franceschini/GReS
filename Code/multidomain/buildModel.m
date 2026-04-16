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

default = struct('Domain',struct(),...
                 'Interface',struct());

params = readInput(default,fileName);
domainInput = params.Domain;
interfInput = params.Interface;

if numel(fieldnames(domainInput))==0
  nD = 0;
else
nD = numel(domainInput);
end

domains = [];

for i = 1:nD
  domains = [domains; defineDomain(domainInput(i))];
end

varargout{1} = domains;

if nargout > 1
  if numel(fieldnames(interfInput))>0
    interfaces = InterfaceSolver.addInterfaces(domains,interfInput);
    varargout{2} = interfaces;
  end

end

end


function domain = defineDomain(input)

input = readInput(input);

topology = Mesh.create(input.Geometry);

if isfield(input,"Materials")
  mat = Materials(input.Materials);
else
  mat = [];
end

gNPoints = 0;

if isfield(input,"Gauss")
  gNPoints = input.Gauss.nGP;
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
  bound = Boundaries(grid,input.BoundaryConditions);
else
  bound = [];
end


domain = Discretizer('grid',grid,...
                     'materials',mat,...
                     'boundaries',bound);

domain.addPhysicsSolvers(input.Solver);


end

