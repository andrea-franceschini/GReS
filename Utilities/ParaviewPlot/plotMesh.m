function plotMesh(mesh, foldName, names)
%PLOTMESH Export one grid or a cell array of grids as a ParaView collection.
%   plotMesh({matrix, fracture}, 'mesh', ["Matrix", "Fracture"]) writes
%   mesh.pvd with one named block per grid. The names are optional.
arguments
  mesh
  foldName (1,1) string
  names string = strings(0,1)
end

if ~iscell(mesh)
  mesh = {mesh};
end
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
