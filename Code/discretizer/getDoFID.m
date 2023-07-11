function dof = getDoFID(mesh,el)
  top = mesh.cells(el,1:mesh.cellNumVerts(el));
  dof = repelem(mesh.nDim*top,mesh.nDim);
  dof = dof + repmat(-(mesh.nDim-1):0,1,mesh.cellNumVerts(el));
end
