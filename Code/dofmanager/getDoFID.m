function dof = getDoFID(mesh,el)
% find multicomponent dofs  in the single physic solassociated
  top = mesh.cells(el,1:mesh.cellNumVerts(el));
  dof = repelem(mesh.nDim*top,mesh.nDim);
  dof = dof + repmat(-(mesh.nDim-1):0,1,mesh.cellNumVerts(el));
end