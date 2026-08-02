function varargout = FEMassembly(solver,localAssembler,varargin)
%FEMASSEMBLY Assemble one or more FEM matrices in bounded-size chunks.
%
% Each call to localAssembler fills the supplied assembler objects.  Their
% sparse outputs are collected in a cell array and added to the contributions
% produced by previous chunks.  One output is returned for each assembler.

default = struct('chunkSize',1e4);
options = readInput(default,varargin{:});
chunkSize = options.chunkSize;

% Number of outputs requested by the caller.
nMat = nargout;

domain = solver.domain;
dofm = domain.dofm;
cells = domain.grid.cells;
subCells = dofm.getFieldCells(solver.getField());
hasContribution = false;

for cTag = 1:cells.nTag

  subRegionCells = subCells(cells.tag(subCells) == cTag);

  for vtkId = cells.vtkTypes
    subCellsLoc = subRegionCells(cells.VTKType(subRegionCells) == vtkId);
    if isempty(subCellsLoc)
      continue
    end

    elem = FiniteElementType.create(vtkId,solver.grid,solver.getGaussOrder);

    for first = 1:chunkSize:numel(subCellsLoc)
      cellList = subCellsLoc(first:min(first + chunkSize - 1,numel(subCellsLoc)));

      % Collect the generic outputs of localAssembler in a cell array.
      matLoc = cell(nMat,1);
      [matLoc{:}] = localAssembler(cTag,cellList,elem);

      if hasContribution
        for m = 1:numel(matLoc)
          mat{m} = mat{m} + matLoc{m};
        end
      else
        % The dimensions and sparse type are known only after the first call.
        for m = 1:numel(matLoc)
          mat = cell(nMat,1);
          mat{m} = matLoc{m};
          hasContribution = true;
        end
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
