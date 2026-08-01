function varargout = FEMassembly(solver,assemblers,localAssembler,varargin)
%FEMASSEMBLY Assemble one or more FEM matrices in bounded-size chunks.
%
% Each call to localAssembler fills the supplied assembler objects.  Their
% sparse outputs are collected in a cell array and added to the contributions
% produced by previous chunks.  One output is returned for each assembler.

default = struct('chunkSize',1e3);
options = readInput(default,varargin{:});
chunkSize = options.chunkSize;

validateattributes(chunkSize,{'numeric'}, ...
  {'scalar','integer','positive'},mfilename,'chunkSize');

domain = solver.domain;
dofm = domain.dofm;
cells = domain.grid.cells;
subCells = dofm.getFieldCells(solver.fieldId);
nMat = numel(assemblers);
mat = cell(nMat,1);
hasContribution = false;

for cTag = 1:cells.nTag
  constLaw = domain.materials.getConstitutiveLaw(cTag);
  subRegionCells = subCells(cells.tag(subCells) == cTag);

  for vtkId = cells.vtkTypes
    subCellsLoc = subRegionCells(cells.VTKType(subRegionCells) == vtkId);
    if isempty(subCellsLoc)
      continue
    end

    elem = FiniteElementType.create(vtkId,solver.grid,solver.gaussOrder);

    for first = 1:chunkSize:numel(subCellsLoc)
      cellList = subCellsLoc(first:min(first + chunkSize - 1,numel(subCellsLoc)));
      nElems = numel(cellList);

      for m = 1:nMat
        assemblers(m).preallocate(nElems);
      end

      localAssembler(assemblers,cellList,elem,constLaw);
      contribution = arrayfun(@(a) a.sparseAssembly(),assemblers, ...
        'UniformOutput',false);
      contribution = contribution(:);

      if hasContribution
        for m = 1:nMat
          mat{m} = mat{m} + contribution{m};
        end
      else
        % The dimensions and sparse type are known only after the first call.
        mat = contribution;
        hasContribution = true;
      end
    end
  end
end

% Preserve the advertised output shape even when there are no active cells.
if ~hasContribution
  for m = 1:nMat
    mat{m} = sparse(assemblers(m).Nrows,assemblers(m).Ncols);
  end
end

varargout = mat;
end
