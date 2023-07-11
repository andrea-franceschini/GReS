classdef Mesh < handle
  % MESH General mesh class, can be subclassed for specialized types

  properties (SetAccess = public, GetAccess = public)

    % GENERAL INFO:
    % Mesh dimension
    nDim = 0
    % Total number of mesh nodes
    nNodes = 0
    % Total number of mesh 3D elements
    nCells = 0
    % Total number of mesh 2D elements
    nSurfaces = 0
    
    % Nodes' coordinates
    coordinates
    % Cell to node mapping:
    % 3D elements' nodes sequences 
    cells
    % 3D elements' tag (region)
    cellTag
    % Number of elements' tag
    nCellTag
    % Number of nodes for each 3D element
    cellNumVerts
    
    % Surface to node mapping:
    % 2D elements' nodes sequences
    surfaces
    % 2D elements' tag (region)
    surfaceTag
    % Number of surfaces' tag
    nSurfaceTag
    % Number of nodes of the 2D element
    surfaceNumVerts

    % Regions
    cellRegions;
    surfaceRegions;

    % 3D element VTK type tag
    cellVTKType;
    % 2D element VTK type tag
    surfaceVTKType;
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

    % Function to read mesh data and store this data inside the properties 
    % of the object created with the class "Mesh"
    function importGMSHmesh(obj, fileName)
      % Reading input file 
      [obj.coordinates, elems, regions] = mxImportGMSHmesh(fileName);
      elems = double(elems);
      
      % STORING DATA INSIDE OBJECT'S PROPERTIES
      % 3D ELEMENT DATA
      obj.nDim = size(obj.coordinates,2);
      obj.nNodes = size(obj.coordinates,1);
      % cellsID = 3D element tag for readGMSHmesh.cpp
      cellsID = [4, 5, 6, 7];
      ID = ismember(elems(:,1), cellsID);
      obj.cellNumVerts = elems(ID,3);
      nVerts = max(obj.cellNumVerts);
      obj.cells = elems(ID,4:nVerts+3);
      obj.cellVTKType = obj.typeMapping(elems(ID,1));
      obj.cellTag = elems(ID,2);
      obj.nCells = length(obj.cellTag);
      %
      % Check for unsupported elements
      if any(~ismember(obj.cellVTKType,[10 12]))
        error(['There are unsupported elements in the mesh.\n', ...
          'Supported elements are: - 4-node tetrahedra (VTKType = %d)\n', ...
          '                        - 8-node hexahedra  (VTKType = %d)'],10,12);
      end
      %
      % REGIONS DATA FOR 3D ELEMENT
      nRegions = length(regions);
      dims = zeros(nRegions,1);
      for i = 1 : nRegions
        dims(i) = regions(i).dim;
      end
      ID = find(dims == 3);
      for i = 1 : length(ID)
        obj.cellRegions = setfield(obj.cellRegions, regions(ID(i)).name, regions(ID(i)).ID);
      end

      % 2D ELEMENT DATA
      % cellsID = 2D surface tag for readGMSHmesh.cpp
      cellsID = [2, 3];
      ID = ismember(elems(:,1), cellsID);
      obj.surfaceNumVerts = elems(ID,3);
      nVerts = max(obj.surfaceNumVerts);
      obj.surfaces = elems(ID,4:nVerts+3);
      obj.surfaceVTKType = obj.typeMapping(elems(ID,1));
      obj.surfaceTag = elems(ID,2);
      obj.nSurfaces = length(obj.surfaceTag);
      %
      % Check for unsupported elements
      if any(~ismember(obj.surfaceVTKType,[5 9]))
        error(['There are unsupported surfaces in the mesh.\n', ...
          'Supported surfaces are: - 3-node triangles       (VTKType = %d)\n', ...
          '                        - 4-node quadrilaterals  (VTKType = %d)'],5,9);
      end
      %
      % REGIONS DATA FOR 2D ELEMENT
      dims = zeros(nRegions,1);
      for i = 1 : nRegions
        dims(i) = regions(i).dim;
      end
      ID = find(dims == 2);
      for i = 1 : length(ID)
        obj.surfaceRegions = setfield(obj.surfaceRegions, regions(ID(i)).name, regions(ID(i)).ID);
      end
      %
      obj.nCellTag = max(obj.cellTag);
      obj.nSurfaceTag = max(obj.surfaceTag);
    end
    
    % Function to call functions that calculate cells' and surfaces' centroids
%     function finalize(obj) 
%       computeCellCentroid(obj);
%       computeSurfaceCentroid(obj);
%     end

    % Function for 3D element centroid calculation
%     function computeCellCentroid(obj)
%       obj.cellCentroid = zeros(obj.nCells,3);
%       for i = 1 : obj.nCells
%         verts = obj.cells(i,1:obj.cellNumVerts(i));
%         obj.cellCentroid(i,:) = sum(obj.coordinates(verts,:),1) / obj.cellNumVerts(i);
%       end
%     end
    %
    
    % Function for 2D element centroid calculation
%     function computeSurfaceCentroid(obj)
%       obj.surfaceCentroid = sparse(repelem(1:obj.nSurfaces,obj.surfaceNumVerts), ...
%           nonzeros((obj.surfaces)'),repelem((obj.surfaceNumVerts).^(-1),obj.surfaceNumVerts),obj.nSurfaces,obj.nNodes) ...
%           * obj.coordinates;
%     end
    %
%     function computeSurfaceCentroid(obj)
%       obj.surfaceCentroid = zeros(obj.nSurfaces,3);
%       for i = 1 : obj.nSurfaces
%         verts = obj.surfaces(i,1:obj.surfaceNumVerts(i));
%         obj.surfaceCentroid(i,:) = sum(obj.coordinates(verts,:),1) / obj.surfaceNumVerts(i);
%       end
%     end

    % Function for call 3D elements' region based on their cellTag
    function ID = findCellsOfRegion(obj, region)
      if (isfield(obj.cellRegions, region))
        val = getfield(obj.cellRegions, region);
      else
        error('Invalid region identifier');
      end
      ID = obj.cellTag == val;
    end

    % Function for call 2D elements' region based on their cellTag
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
