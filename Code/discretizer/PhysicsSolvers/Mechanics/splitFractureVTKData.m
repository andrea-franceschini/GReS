function [bulk,bulkData,fractures,fractureData] = splitFractureVTKData(model,jump,trace)
% Export rock-side displacements and facetwise jump fields. Use a separate
% point per facet vertex so a crossing can carry several different jumps.
dim = size(model.points,2);
bulkType = [5,10];
facetType = [3,5];
bulk = struct('points',padXYZ(model.points(model.traceVertex,:)), ...
  'cells',model.simplexTraces,'vtkType',bulkType(dim-1));
bulkData = struct('name','displacements','data',padXYZ(trace));

nFacets = size(model.facets,1);
vertices = reshape(model.facets.',[],1);
fractures = struct('points',padXYZ(model.points(vertices,:)), ...
  'cells',reshape(1:numel(vertices),dim,nFacets).','vtkType',facetType(dim-1));
jump = reshape(permute(jump,[1,3,2]),[],dim);
normal = repelem(model.normal,dim,1);
opening = sum(jump.*normal,2);
slip = jump-opening.*normal;
fractureData = struct('name',{'fractureJump';'normalOpening';'tangentialJump';'fracturePressure'}, ...
  'data',{padXYZ(jump);opening;padXYZ(slip);repelem(model.pressure,dim,1)});
end

function xyz = padXYZ(values)
xyz = [values,zeros(size(values,1),3-size(values,2))];
end
