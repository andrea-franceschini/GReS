function grdecl = refineLayers(grdecl,zMod,nRefs)

% input is a corner point grid as obtained from mrst

% the output grids refine the z later (from below) specified in zMode.
% refRatio: the number of subdivision for the input layer

% identify starting index of zcorns for the selected layer.
assert(zMod < grdecl.cartDims(3),'Layer to refine exceed number of available vertical layers');

nPointsLayer = 4*grdecl.cartDims(1)*grdecl.cartDims(2);

if zMod == 1
  k = 1;
else
  k = 2*zMod-1;
end

ptBot = nPointsLayer*(k-1);

zBot = grdecl.ZCORN(ptBot+1:ptBot+nPointsLayer);
zTop = grdecl.ZCORN(ptBot+nPointsLayer+1:ptBot+2*nPointsLayer);

dz = zTop-zBot;
zNew = zeros(numel(dz),2*nRefs);

for i = 1:length(zBot)
  zVals = linspace(zBot(i),zTop(i),nRefs+1);
  zNew(i,:) = [zVals(1), repelem(zVals(2:end-1),2), zVals(end)];
end

grdecl.ZCORN = [grdecl.ZCORN(1:ptBot);
                zNew(:);
                grdecl.ZCORN(ptBot+2*nPointsLayer+1:end)];


% update rest of grid structure
grdecl.cartDims(3) = grdecl.cartDims(3) + nRefs -1;
numbNewCells = grdecl.cartDims(1)*grdecl.cartDims(2)*(nRefs-1);
grdecl.ACTNUM = [grdecl.ACTNUM; ones(numbNewCells,1)];

% check correct number of zcorn
assert(length(grdecl.ZCORN)==8*prod(grdecl.cartDims))


% add refRatio layers with zcorn dividing the bounding coordinates

% update the number of active cells and the number of layers



end

