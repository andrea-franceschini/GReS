function [coordinates, tetrahedra, fractures] = meshFracturedBox(bounds, polygons, ...
    matrixSize, fractureSize, transitionFactor)
%MESHFRACTUREDBOX Return a conforming tetrahedral mesh and per-fracture triangles.
for k = 1:numel(polygons)
  assert(all(all(polygons{k} >= bounds(:,1).' & polygons{k} <= bounds(:,2).')), ...
    'GReS:GmshInput', 'Fracture vertices must lie in the box.');
end

gresgmsh.call('initialize');
session = onCleanup(@() gresgmsh.call('finalize')); %#ok<NASGU>
gresgmsh.call('option.setNumber', 'General.Terminal', 0);
gresgmsh.call('model.add', 'fracturedBox');
origin = num2cell(bounds(:,1));
lengths = num2cell(diff(bounds,1,2));
box = gresgmsh.call('model.occ.addBox', origin{:}, lengths{:}, -1);

fractureSurfaces = cell(size(polygons));
if ~isempty(polygons)
  tools = zeros(2,numel(polygons),'int32');
  for k = 1:numel(polygons)
    tools(:,k) = [2; addPolygon(polygons{k})];
  end
  [~, mapping] = gresgmsh.call('model.occ.fragment', int32([3 box]), tools(:), -1, 1, 1);
  assert(numel(mapping) == numel(polygons) + 1, 'GReS:GmshFragmentMap', ...
    'Unexpected OCC fragment map size.');
  for k = 1:numel(polygons)
    entities = reshape(mapping{k+1},2,[]);
    fractureSurfaces{k} = entities(2,entities(1,:) == 2);
  end
end
gresgmsh.call('model.occ.synchronize');
gresgmsh.call('model.mesh.generateGraded', polygons, matrixSize, fractureSize, transitionFactor);

[nodeTags, xyz] = gresgmsh.call('model.mesh.getNodes');
assert(~isempty(nodeTags), 'GReS:GmshMesh', 'Gmsh returned no nodes.');
coordinates = reshape(xyz,3,[]).';

tet = gresgmsh.call('model.mesh.getElementType', 'Tetrahedron', 1, 0);
[~, tags] = gresgmsh.call('model.mesh.getElementsByType', tet, -1);
assert(~isempty(tags), 'GReS:GmshMesh', 'Gmsh returned no tetrahedra.');
tetrahedra = connectivity(tags, nodeTags, 4);

fractures = cell(numel(polygons),1);
if ~isempty(polygons)
  tri = gresgmsh.call('model.mesh.getElementType', 'Triangle', 1, 0);
  for k = 1:numel(polygons)
    surfaces = double(fractureSurfaces{k});
    blocks = cell(size(surfaces));
    for j = 1:numel(surfaces)
      [~, blocks{j}] = gresgmsh.call('model.mesh.getElementsByType', tri, surfaces(j));
    end
    tags = [blocks{:}];
    assert(~isempty(tags), 'GReS:GmshMesh', 'Gmsh returned no triangles for fracture %d.',k);
    fractures{k} = connectivity(tags, nodeTags, 3);
  end
end
end

function surface = addPolygon(polygon)
n = size(polygon,1);
points = zeros(n,1);
lines = zeros(n,1,'int32');
for k = 1:n
  points(k) = gresgmsh.call('model.occ.addPoint', polygon(k,1), polygon(k,2), polygon(k,3), 0, -1);
end
for k = 1:n
  lines(k) = gresgmsh.call('model.occ.addLine', points(k), points(mod(k,n)+1), -1);
end
loop = gresgmsh.call('model.occ.addCurveLoop', lines, -1);
surface = gresgmsh.call('model.occ.addPlaneSurface', int32(loop), -1);
end

function indices = connectivity(tags, nodeTags, width)
[found, indices] = ismember(tags,nodeTags);
assert(all(found), 'GReS:GmshMesh', 'Element connectivity references an unknown node.');
indices = reshape(indices,width,[]).';
end
