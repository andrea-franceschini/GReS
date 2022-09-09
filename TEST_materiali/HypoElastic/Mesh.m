classdef Mesh < handle
  % MESH General mesh class, can be subclassed for specialized types

  properties (SetAccess = public, GetAccess = public)

    % GENERAL INFO:
    % Mesh dimension
    nDim = 0; 
    % Total number of mesh nodes
    nNodes = 0;
    % Total number of mesh 3D elements
    nCells = 0;
    % Total number of mesh 2D elements
    nSurfaces = 0;

    % Nodal coordinates
    coordinates = [];

    % Cell to node mapping:
    % Nodes of i-element
    cells = [];
    %Material tag
    cellTag = [];
    %Number of nodes of the 3D element
    cellNumVerts = [];
    
    %cellToNode = logical(sparse(0, 0));
    %nodeToCell = logical(sparse(0, 0));

    % Surface to node mapping:
    % Nodes of i-surface
    surfaces = [];
    % Surface tag (position)
    surfaceTag = [];
    % Number of nodes of the 2D element
    surfaceNumVerts = [];

    % Regions
    cellRegions = [];
    surfaceRegions = [];

    % Coordinates of 3D element centroid
    cellCentroid = [];
    % Coordinates of 2D element centroid
    surfaceCentroid = [];

    % Tag for 3D element type
    cellVTKType = [];
    % Tag for 2D element type
    surfaceVTKType = [];
    meshType = 'Unstructured'

  end

  properties (Access = private)

    %  1 VTK_VERTEX
    %  2 VTK_POLY_VERTEX
    %  3 VTK_LINE
    %  4 VTK_POLY_LINE
    %  5 VTK_TRIANGLE
    %  6 VTK_TRIANGLE_STRIP
    %  7 VTK_POLYGON
    %  8 VTK_PIXEL
    %  9 VTK_QUAD
    % 10 VTK_TETRA
    % 11 VTK_VOXEL
    % 12 VTK_HEXAHEDRON
    % 13 VTK_WEDGE
    % 14 VTK_PYRAMID
    typeMapping = [3; 5; 9; 10; 12; 13; 14];

  end

  methods (Access = public)

    function importGMSHmesh(obj, fileName)
      % Reading input file, creation of 3 matrices
      [obj.coordinates, elems, regions] = mxImportGMSHmesh(fileName);
      % Conversion to double precision
      elems = double(elems);
      % Assignment of prop nDim(scalar) = number of columns of obj.coordinates
      obj.nDim = size(obj.coordinates,2);
      % Assignment of prop nNodes(scalar) = number of rows of obj.coordinates
      obj.nNodes = size(obj.coordinates,1);
      % Vector cellsID: 3D element tag for readGMSHmesh.cpp
      cellsID = [4, 5, 6, 7];
      % Vector ID: for 3D elements' position in matrix "elems"
      ID = ismember(elems(:,1), cellsID);
      % Assignment of prop cellNumVerts(vector,size = nCells) = number of the nodes of the 3D elements
      obj.cellNumVerts = elems(ID,3);
      % Assignment of nVerts (scalar) = max element of obj.cellNumVerts
      nVerts = max(obj.cellNumVerts);
      % Assignment of prop cells (vector) = columns 4 to nVerts+3 of "elems"
      obj.cells = elems(ID,4:nVerts+3);
      % Assignment of prop cellVTKtype(vector) = first column of "elems" 
      obj.cellVTKType = obj.typeMapping(elems(ID,1));
      % Assignment of prop cellTag (vector) = second column of "elems"
      obj.cellTag = elems(ID,2);
      % Assignment of prop nCells (scalar) = size of cellTag
      obj.nCells = length(obj.cellTag);

      % Creation of index nRegions
      nRegions = length(regions);
      % Space allocation for vector "dims"
      dims = zeros(nRegions,1);
      for i = 1 : nRegions
        dims(i) = regions(i).dim;
      end
      % Creation of vector with linear indices of =3 entries in dims
      ID = find(dims == 3);
      for i = 1 : length(ID)
        obj.cellRegions = setfield(obj.cellRegions, regions(ID(i)).name, regions(ID(i)).ID);
      end

      % Same as cells' prop, but with surfaces
      %vector cellsID: 2D element tag for readGMSHmesh.cpp
      cellsID = [2, 3];
      %vector ID: for 2Delements' position in "elems"
      ID = ismember(elems(:,1), cellsID);
      obj.surfaceNumVerts = elems(ID,3);
      nVerts = max(obj.surfaceNumVerts);
      obj.surfaces = elems(ID,4:nVerts+3);
      obj.surfaceVTKType = obj.typeMapping(elems(ID,1));
      obj.surfaceTag = elems(ID,2);
      obj.nSurfaces = length(obj.surfaceTag);

      dims = zeros(nRegions,1);
      for i = 1 : nRegions
        dims(i) = regions(i).dim;
      end
      % Creation of vector with linear indices of =2 entries in dims
      ID = find(dims == 2);
      for i = 1 : length(ID)
        obj.surfaceRegions = setfield(obj.surfaceRegions, regions(ID(i)).name, regions(ID(i)).ID);
      end
    end
    
    function finalize(obj) 
      computeCellCentroid(obj);
      computeSurfaceCentroid(obj);
    end

    % Function for 3D element centroid calculation
    function computeCellCentroid(obj)
      obj.cellCentroid = zeros(obj.nCells,3);
      for i = 1 : obj.nCells
        verts = obj.cells(i,1:obj.cellNumVerts(i));
        obj.cellCentroid(i,:) = sum(obj.coordinates(verts,:),1) / obj.cellNumVerts(i);
      end
    end
    
    % Function for 2D element centroid calculation
    function computeSurfaceCentroid(obj)
      obj.surfaceCentroid = zeros(obj.nSurfaces,3);
      for i = 1 : obj.nSurfaces
        verts = obj.surfaces(i,1:obj.surfaceNumVerts(i));
        obj.surfaceCentroid(i,:) = sum(obj.coordinates(verts,:),1) / obj.surfaceNumVerts(i);
      end
    end

    function ID = findCellsOfRegion(obj, region)
      if (isfield(obj.cellRegions, region))
        val = getfield(obj.cellRegions, region);
      else
        error('Invalid region identifier');
      end
      ID = obj.cellTag == val;
    end

    function ID = findSurfacesOfRegion(obj, region)
      if (isfield(obj.surfaceRegions, region))
        val = getfield(obj.surfaceRegions, region);
      else
        error('Invalid region identifier');
      end
      ID = obj.surfaceTag == val;
    end

  end

end
