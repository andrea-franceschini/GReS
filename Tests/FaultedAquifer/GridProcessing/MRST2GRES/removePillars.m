function grdecl = removePillars(grdecl, dir, cellsToRemove)

% input is a corner point grid as obtained from mrst
grdeclOriginal  = grdecl;

idPillars = 1:prod(grdecl.cartDims(1:2)+1);
idPillars = reshape(idPillars,grdecl.cartDims(1:2)+1);
if dir==1
  idPillars(cellsToRemove,:) = [];
else
  idPillars(:,cellsToRemove) = [];
end

idPillars = idPillars(:);
idPillars = repelem((idPillars-1)*6, 6) + repmat((1:6)', numel(idPillars), 1);


idCells = 1:8*prod(grdecl.cartDims);
idCells = reshape(idCells,2*grdecl.cartDims);

idC = 2*cellsToRemove(1)-1:2*cellsToRemove(end);

% remove cells
if dir==1
  idCells(idC,:,:) = [];
else
  idCells(:,idC,:) = [];
end

idCells = idCells(:);


grdecl.COORD = grdecl.COORD(idPillars);
grdecl.ZCORN = grdecl.ZCORN(idCells);

grdecl.cartDims(dir) = grdecl.cartDims(dir) - numel(cellsToRemove);
grdecl.ACTNUM = int32(ones(prod(grdecl.cartDims),1));



end








