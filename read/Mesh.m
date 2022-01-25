classdef Mesh < handle
  % MESH General mesh class, can be subclassed for specialized types

  properties (SetAccess = public, GetAccess = public)

    % general info
    nDim = 0;
    nNodes = 0;
    nCells = 0;
    nSurfaces = 0;

    % nodal coordinates
    coordinates = [];

    % cell to node mapping
    cells = [];
    cellTag = [];
    cellNumVerts = [];
    %cellToNode = logical(sparse(0, 0));
    %nodeToCell = logical(sparse(0, 0));

    % surface to node mapping
    surfaces = [];
    surfaceTag = [];
    surfaceNumVerts = [];

    % Regions
    cellRegions = [];
    surfaceRegions = [];

    cellCentroid = [];
    surfaceCentroid = [];

    cellVTKType = [];
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
      [obj.coordinates, elems, regions] = mxImportGMSHmesh(fileName);
      elems = double(elems);
      obj.nDim = size(obj.coordinates,2);
      obj.nNodes = size(obj.coordinates,1);
      cellsID = [4, 5, 6, 7];
      ID = ismember(elems(:,1), cellsID);
      obj.cellNumVerts = elems(ID,3);
      nVerts = max(obj.cellNumVerts);
      obj.cells = elems(ID,4:nVerts+3);
      obj.cellVTKType = obj.typeMapping(elems(ID,1));
      obj.cellTag = elems(ID,2);
      obj.nCells = length(obj.cellTag);

      nRegions = length(regions);
      dims = zeros(nRegions,1);
      for i = 1 : nRegions
        dims(i) = regions(i).dim;
      end
      ID = find(dims == 3);
      for i = 1 : length(ID)
        obj.cellRegions = setfield(obj.cellRegions, regions(ID(i)).name, regions(ID(i)).ID);
      end

      cellsID = [2, 3];
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
      ID = find(dims == 2);
      for i = 1 : length(ID)
        obj.surfaceRegions = setfield(obj.surfaceRegions, regions(ID(i)).name, regions(ID(i)).ID);
      end
    end

    function finalize(obj)
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
