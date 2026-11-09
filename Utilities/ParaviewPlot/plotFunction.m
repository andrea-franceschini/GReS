function plotFunction(mesh, foldName, time, funct, location)
%PLOTFUNCTION Export a nodal or elemental field as a ParaView collection.
%   plotFunction(mesh, 'solution', time, values, 'elem') writes cell data.
%   The default location is 'node'.
arguments
  mesh
  foldName (1,1) string
  time (1,1) double {mustBeReal,mustBeFinite}
  funct
  location (1,1) string = "node"
end

location = lower(location);
assert(any(location == ["node","elem"]), 'GReS:PlotFunction', ...
  'Field location must be node or elem.');
[~, name] = fileparts(foldName);
writePlotCollection({mesh},foldName,time,{funct},location,name);
end
