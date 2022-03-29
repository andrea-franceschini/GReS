classdef Mesh < handle
  % MESH General mesh class, can be subclassed for specialized types

  properties (SetAccess = public, GetAccess = public)

    % general info
    nDim = 0; %Mesh dimension
    nNodes = 0;%Total number of mesh nodes
    nCells = 0;%Total number of mesh 3D elements
    nSurfaces = 0;%Total number of mesh 2D elements

    % nodal coordinates
    coordinates = [];

    % cell to node mapping
    cells = [];%nodes of i-element
    cellTag = [];%Material tag
    cellNumVerts = [];%Number of nodes of the   3D element
    %cellToNode = logical(sparse(0, 0));
    %nodeToCell = logical(sparse(0, 0));

    % surface to node mapping
    surfaces = [];%nodes of i-surface
    surfaceTag = [];%Surface tag (position)
    surfaceNumVerts = [];%Number of nodes of the 2D element

    % Regions
    cellRegions = [];
    surfaceRegions = [];

    cellCentroid = [];%coordinates of 3D element centroid
    surfaceCentroid = [];%coordinates of 2D element centroid

    cellVTKType = [];%Tag for 3Delement type
    surfaceVTKType = [];%Tag for 2D element type
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
      [obj.coordinates, elems, regions] = mxImportGMSHmesh(fileName);%reading input file, creation of 3 matrices
      elems = double(elems);%conversion to double precision
      obj.nDim = size(obj.coordinates,2);%assignment of prop nDim(scalar) = number of columns of obj.coordinates 
      obj.nNodes = size(obj.coordinates,1);%assignment of prop nNodes(scalar) = number of rows of obj.coordinates
      cellsID = [4, 5, 6, 7];%vector with 3D element tag for readGMSHmesh.cpp
      ID = ismember(elems(:,1), cellsID);%vector for 3Delements' position in matrix "elems" 
      obj.cellNumVerts = elems(ID,3);%assignment of prop cellNumVerts(vector,size = nCells) = number of the nodes of the 3D elements 
      nVerts = max(obj.cellNumVerts);%assignment of nVerts (scalar) = max element of obj.cellNumVerts
      obj.cells = elems(ID,4:nVerts+3);%assignment of prop cells (vector) = columns 4 to nVerts+3 of "elems"
      obj.cellVTKType = obj.typeMapping(elems(ID,1));%assignment of prop cellVTKtype(vector) = first column of "elems" 
      obj.cellTag = elems(ID,2);%assignment of prop cellTag (vector) = second column of "elems"
      obj.nCells = length(obj.cellTag);%assignment of prop nCells (scalar) = size of cellTag

      nRegions = length(regions);%creation of index nRegions
      dims = zeros(nRegions,1);%dims allocation vector
      for i = 1 : nRegions
        dims(i) = regions(i).dim;
      end
      ID = find(dims == 3);%creation of vector with linear indices of =3 entries in dims
      for i = 1 : length(ID)
        obj.cellRegions = setfield(obj.cellRegions, regions(ID(i)).name, regions(ID(i)).ID);%cellRegions.name(i) = ID(i)?
      end

      %same as cells' prop, but with surfaces
      cellsID = [2, 3];%vector with 2D element tag for readGMSHmesh.cpp
      ID = ismember(elems(:,1), cellsID);%vector for 2Delements' position in "elems"
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
      ID = find(dims == 2);%creation of vector with linear indices of =2 entries in dims
      for i = 1 : length(ID)
        obj.surfaceRegions = setfield(obj.surfaceRegions, regions(ID(i)).name, regions(ID(i)).ID);
      end
    end

    function finalize(obj) %function for centroids  calculation
      computeCellCentroid(obj);
      computeSurfaceCentroid(obj);
    end

    function computeCellCentroid(obj)
      obj.cellCentroid = zeros(obj.nCells,3);
      for i = 1 : obj.nCells
        verts = obj.cells(i,1:obj.cellNumVerts(i));
        obj.cellCentroid(i,:) = sum(obj.coordinates(verts,:),1) / obj.cellNumVerts(i);
      end
    end

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
