function plotMesh(mesh, foldName, names)
%PLOTMESH Export one grid or a cell array of grids as a ParaView collection.
%   plotMesh({matrix, fracture}, 'mesh', ["Matrix", "Fracture"]) writes
%   mesh.pvd with one named block per grid. The names are optional.
%   A fracture cell array from createFracturedGrid can be passed as a block list.
narginchk(2,3);
foldName = string(foldName);
assert(isscalar(foldName), 'GReS:PlotMesh', 'Provide one output folder name.');
if nargin < 3
  names = strings(0,1);
else
  names = string(names);
end

if ~iscell(mesh)
  mesh = {mesh};
end
grids = {};
for k = 1:numel(mesh)
  if iscell(mesh{k})
    grids = [grids, mesh{k}(:).']; %#ok<AGROW>
  else
    grids{end+1} = mesh{k}; %#ok<AGROW>
  end
end
mesh = grids;
assert(~isempty(mesh), 'GReS:PlotMesh', 'At least one grid is required.');
if isempty(names)
  if isscalar(mesh)
    [~, names] = fileparts(foldName);
  else
    names = "Grid " + (1:numel(mesh));
  end
end
assert(numel(names) == numel(mesh) && all(strlength(names(:)) > 0), ...
  'GReS:PlotMesh', 'Provide one nonempty block name per grid.');

values = cellfun(@(grid) zeros(grid.nNodes,1),mesh,'UniformOutput',false);
writePlotCollection(mesh,foldName,0,values,"node",names);
end
