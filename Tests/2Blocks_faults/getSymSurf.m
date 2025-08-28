function list = getSymSurf(interf)
switch interf.multiplierType
  case 'P0'
    % get list of surfaces along the symmetry axis of the fault
    % it is needed for the print utilities of the 2 blocks fault model
    msh = interf.mesh.msh(2);
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
  otherwise
    % get list of surfaces along the symmetry axis of the fault
    % it is needed for the print utilities of the 2 blocks fault model
    msh = interf.mesh.msh(2);
    % get surfaces on the symmetry axis
    tol = 1e-3; % Set a small tolerance
    id = abs(msh.coordinates(:,2)-5)<tol;
    % get cell list exploiting mesh regularity
    list = find(id);
    % sort depending on the z- coordinate
    [~,idsort] = sort(msh.coordinates(list,3));
    list = list(idsort);
end
end

