function FEMassembly(solver,localAssembler)

% localAssembler is an anonymous function requiring:
% a list of cells
% the constitutiveLaw object
% the finite element object

% routine that loops in a domain and call a local assembler on a subset of
% cells having same element type and region tag

domain = solver.domain;

dofm = domain.dofm;
cells = domain.grid.cells;

subCells = dofm.getFieldCells(obj.fieldId);

for cTag = 1:cells.nTag

  % extract the constitutive law
  constLaw = obj.domain.materials.getConstitutiveLaw(cTag);

  % extract cells belonging to subregion
  subRegionCells = subCells(cells.tag(subCells) == cTag);


  for vtkId = cells.vtkTypes

    % extract cells of the subregion of homogeneous vtk type
    cellId = cells.VTKType(subRegionCells) == vtkId;
    subCellsLoc = subRegionCells(cellId);

    cellList = find(cellId);

    elem = FiniteElementType.create(vtkId,obj.grid,obj.gaussOrder);

    localAssembler()



  end
end

end

