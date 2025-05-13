function list = getSymSurf(mG)
% get list of surfaces along the symmetry axis of the fault
% it is needed for the print utilities of the 2 blocks fault model
msh = mG.interfaces.mortar.intSlave;
% get surfaces on the symmetry axis
tolerance = 1e-3; % Set a small tolerance
x_coords = unique(round(msh.coordinates(:,2)/tolerance)*tolerance);
y_coords = unique(round(msh.coordinates(:,3)/tolerance)*tolerance); 
nxCells = length(x_coords)-1;
nyCells = length(y_coords)-1;
% get cell list exploiting mesh regularity
nx = round(nxCells/2);
l1 = (nx-1)*nyCells+1; 
l2 = nx*nyCells;
list = (l1:l2)';
end

