function varargout = FEMassembly(solver,localAssembler,varargin)
%FEMASSEMBLY 
% Efficient Kernel for performing finite element computations

default = struct('chunkSize',1e4);
options = readInput(default,varargin{:});
chunkSize = options.chunkSize;

nMat = nargout;
mat = cell(nMat,1);
hasContribution = false;

domain = solver.domain;
dofm = domain.dofm;
cells = domain.grid.cells;
subCells = dofm.getFieldCells(solver.getField());

for cTag = 1:cells.nTag

  subRegionCells = subCells(cells.tag(subCells) == cTag);

  for vtkId = cells.vtkTypes

    subCellsLoc = ...
      subRegionCells(cells.VTKType(subRegionCells) == vtkId);

    if isempty(subCellsLoc)
      continue
    end

    elem = FiniteElementType.create( ...
      vtkId,solver.grid,solver.getGaussOrder);

    for first = 1:chunkSize:numel(subCellsLoc)

      cellList = subCellsLoc( ...
        first:min(first + chunkSize - 1,numel(subCellsLoc)));

      if nMat == 0
        % Execute localAssembler only for its side effects.
        localAssembler(cellList,elem);
        continue
      end

      matLoc = cell(nMat,1);
      [matLoc{:}] = localAssembler(cellList,elem);

      if ~hasContribution
        mat = matLoc;
        hasContribution = true;
      else
        for m = 1:nMat
          mat{m} = mat{m} + matLoc{m};
        end
      end

    end
  end
end

if nMat > 0 && ~hasContribution
  error('FEMassembly:noContributions', ...
    ['No active cells were found, so the dimensions of the output ', ...
     'matrices could not be inferred.']);
end

varargout = mat;

end