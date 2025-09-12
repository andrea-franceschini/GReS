function vals = setTabularParams(fName,mesh)
% read material parameter file
% neirest neighbor interpolation is adopted to ensure that all cells
% have a parameter assigned
constVal = str2num(fName);
if ~isempty(constVal)
   vals = constVal*ones(mesh.nCells,1);
   return
end
fID = Materials.openReadOnlyFile(fName);
nC = [];
while isempty(nC)
   line = readToken(fID,fName);
   if ~isempty(line) && ~strcmp(line(1),'%')
      nC = str2double(line);
   else
      continue
   end
end
assert(nC<=mesh.nCells,['Too many input cells in ' ...
   'file %s'],fName);
availCells = zeros(nC,1);
availVals = zeros(nC,1);
% read cell list
c = 0;
while c < nC
   line = fgetl(fID);
   if isempty(line) || strcmp(line(1),'%')
      continue
   end
   availCells(c+1) = sscanf(line, '%i');
   c = c+1;
end
assert(max(availCells)<mesh.nCells,['Cell index out of available cell in ' ...
   'file %s'],fName);
assert(numel(unique(availCells))==numel(availCells),['Repeated cell index' ...
   ' in file %s'],fName)
% read values list
c = 0;
while c < nC
   line = fgetl(fID);
   if isempty(line) || strcmp(line(1),'%')
      continue
   end
   availVals(c+1) = sscanf(line, '%e');
   c = c+1;
end
assert(feof(fID),'Too man parameters value in %s',fName);

if numel(availVals) == mesh.nCells
   % neirest neighbor interpolation not needed
   vals = availVals;
else
   vals = zeros(mesh.nCells,1);
   vals(availCells) = availVals;
   cID = find(~ismember(1:mesh.nCells,availCells));
   for j = cID
      vals(j) = neirestNeigh(j,mesh,availCells,availVals);
   end
end
end

function val = neirestNeigh(id,mesh,cells,vals)
   % compute distance of cell centroid against available cell centroids;
   d = sqrt(sum((mesh.cellCentroid(cells,:)-mesh.cellCentroid(id,:)).^2,2));
   [~,id] = min(d);
   val = vals(id);
end
