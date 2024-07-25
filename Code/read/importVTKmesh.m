function [coord,elems] = importVTKmesh(fileName)
% Read VTK mesh calling MEX-file and post-processing output structure
% Cell data handling still missing

out = mx_vtkRead(fileName); % call to mex file

% collecting data struct in GReS format
coord = out.points;
nV = sum(out.cells>0,2);
nC = length(out.cellTypes);
elems = zeros(nC,3+max(nV));
elems(:,1) = out.cellTypes;
if isempty(fieldnames(out.cellData)) % cellData not specified
elems(:,2) = 1;
end
elems(:,3) = nV;
elems(:,4:end) = out.cells;
end

