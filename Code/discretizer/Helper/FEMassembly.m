function varargout = FEMassembly(solver,assemblers,localAssembler,chunkSize)


% assemblers is an array of preallocated assembler objects corresponding to
% the output matrices
% localAssembler is an anonymous function requiring:
% a list of cells
% the constitutiveLaw object
% the finite element object

% optional: chunkSize -> the maximum number of cells processed together,
% default is 1e3.

% output: a matrix for each of the input assembler objects.



% routine that loops in a domain and call a local assembler on a subset of
% cells having same element type and region tag

if nargin < 4
  chunkSize = 1e3;
end

domain = solver.domain;

dofm = domain.dofm;
cells = domain.grid.cells;

subCells = dofm.getFieldCells(obj.fieldId);

nMat = numel(assemblers);
mat = cell(nMat,1);

% prepare sparse matrices

for m = 1:nMat


  mat{m} = sparse(assemblers(m).Nrows,assemblers(m).Ncols);



end

for cTag = 1:cells.nTag

  % extract the constitutive law
  constLaw = obj.domain.materials.getConstitutiveLaw(cTag);

  % extract cells belonging to subregion
  subRegionCells = subCells(cells.tag(subCells) == cTag);


  for vtkId = cells.vtkTypes

    % extract cells of the subregion of homogeneous vtk type
    cellId = cells.VTKType(subRegionCells) == vtkId;
    subCellsLoc = subRegionCells(cellId);

    elem = FiniteElementType.create(vtkId,obj.grid,obj.gaussOrder);

    k = 0;

    nChunks = floor(numel(subCellsLoc) / chunkSize);
    lastChunk = rem(numel(subCellsLoc), chunkSize);
    chunks = repmat(chunkSize,nChunks);
    if lastChunk > 0
      chunks(end+1) = lastChunk;
      nChunks = nChunks + 1;
    end


    k = 0;

    for i = 1:nChunks

      nElems = chunks(i);
      cellList = subCellsLoc(k+1:k+nElems);

      % element wise operation on assemblers. call preallocate

      for m = 1:nMat

        assemblers(m).preallocate(nElems);

      end

      % call assembler for local elements
      localAssembler(cellList,elem,constLaw,assemblers)

      for m = 1:nMat


        mat{m} = mat{m} + assemblers(m).sparseAssembly;


      end

      k = k + nElems;



    end


  end
end


varargout = mat;

end

