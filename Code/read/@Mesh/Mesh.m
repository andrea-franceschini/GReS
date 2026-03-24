classdef Mesh < handle


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
    % Total number of mesh 1D elements
    nEdges = 0
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
    % Centroid coordinates of each cell
    cellCentroid
    % Centroid coordinates of each cell
    surfaceCentroid
    % volume of cells
    cellVolume
    % area of surfaces
    surfaceArea
    % Surface to node mapping:
    % 2D elements' nodes sequences
    surfaces
    % 2D elements' tag (region)
    surfaceTag
    % Number of surfaces' tag
    nSurfaceTag
    % Number of nodes of the 2D element
    surfaceNumVerts
    % Edge to node mapping (Required for 2D problems):
    % 1D elements' nodes sequences
    edges
    % 1D elements' tag (region)
    edgeTag
    % Number of nodes of the 2D element
    edgeNumVerts
    % Flag for Cartesian grids
    cartGrid = false;
    % Regions
    cellRegions;
    surfaceRegions;
    % 3D element VTK type tag
    cellVTKType;
    % 2D element VTK type tag
    surfaceVTKType;
    % 1D element VTK type tag
    edgeVTKType
    %
    meshType = 'Unstructured'
    %meshType = 'Undefined'
  end

  properties (Access = private)

    %  1  VTK_VERTEX
    %  2  VTK_POLY_VERTEX
    %  3  VTK_LINE
    %  4  VTK_POLY_LINE
    %  5  VTK_TRIANGLE
    %  6  VTK_TRIANGLE_STRIP
    %  7  VTK_POLYGON
    %  8  VTK_PIXEL
    %  9  VTK_QUAD
    % 10  VTK_TETRA
    % 11  VTK_VOXEL
    % 12  VTK_HEXAHEDRON
    % 13  VTK_WEDGE
    % 14  VTK_PYRAMID
    % 28  VTK_BIQUADRATIC_QUAD
    % 29  VTK_TRIQUADRATIC_HEXAHEDRON

    % Available types in gmsh
    cellVTK = [10, 12, 29];
    surfaceVTK = [5, 9, 28];
    edgeVTK = [3,21];
  end

  methods (Access = public)

    importMesh(obj,fileName)

  
    function [surfMesh,varargout] = getSurfaceMesh(obj, surfTag, varargin)
        % Function to build a separate 2D mesh object based on the
        % surfaceTag of another object
        % there are 3 was to call this method
        % 1) in1 = surfaceTag(s)
        % 2) in1 = surfaceTags in2= logical index of active surfaceTags
        % 3) in1 = logical index of active surfaces

        % optional output: global node indices of surface mesh
        % initialize Mesh object
        if nargin>2
          assert(~isempty(varargin{1}));
          id = find(ismember(obj.surfaceTag,surfTag));
          id = id(varargin{1});
        else
          if ~islogical(surfTag)
            id = ismember(obj.surfaceTag,surfTag);
          else
            assert(numel(surfTag) == obj.nSurfaces,['Logical ' ...
              'index array must have the same size of available surfaces'])
            id = surfTag;
          end
        end
        if sum(id)==0
          surfMesh = [];
          return
        end
        surfMesh = Mesh();
        surfTopol = obj.surfaces(id,:);
        if nargout > 1
          % global node index ordered per column
          varargout{1} = surfTopol;
        end
        % renumber the nodes starting from 1;
        surfTopol = surfTopol(:);
        % ordered list of unique nodes in the topology matrix
        [surfOrd] = unique(surfTopol);
        if surfOrd(1) == 0
          val = 0:numel(surfOrd)-1;
        else
          val = 1:numel(surfOrd);
        end
        mapping = containers.Map(surfOrd, val);
        for i = 1:length(surfTopol)
            surfTopol(i) = mapping(surfTopol(i));
        end
        if surfOrd(1) == 0
          surfOrd = surfOrd(2:end);
        end

        surfMesh.surfaceNumVerts = Mesh.copyField(obj.surfaceNumVerts,id);
        nNmax = max(surfMesh.surfaceNumVerts);
        surfMesh.surfaces = (reshape(surfTopol, [], nNmax));
        surfMesh.coordinates = obj.coordinates(surfOrd,:);
        surfMesh.nNodes = length(surfMesh.coordinates);
        surfMesh.nSurfaces = size(surfMesh.surfaces,1);
        surfMesh.surfaceTag = obj.surfaceTag(id);
        surfMesh.nSurfaceTag = 1;
        surfMesh.surfaceVTKType = Mesh.copyField(obj.surfaceVTKType,id);
        surfMesh.nDim = 3;
        surfMesh.surfaceCentroid = Mesh.copyField(obj.surfaceCentroid,id);
        surfMesh.surfaceArea = Mesh.copyField(obj.surfaceArea,id);
        surfMesh.meshType = obj.meshType;
    end


    function addSurface(obj,id,topol)
       % add a surface to mesh object given the surface topology
       surf = load(topol); % standard topology file
       obj.surfaces = [obj.surfaces; surf];
       if any(obj.surfaceTag==id)
          error('Surface iD already been defined')
       else
          obj.surfaceTag = [obj.surfaceTag; id*ones(size(surf,1),1)];
       end

       if isempty(obj.nSurfaceTag)
          obj.nSurfaceTag = 1;
       else
          obj.nSurfaceTag = obj.nSurfaceTag + 1;
       end
       obj.surfaceNumVerts = [obj.surfaceNumVerts; sum(surf > 0,2)];
       obj.surfaceVTKType(obj.surfaceNumVerts == 3) = 5;
       obj.surfaceVTKType(obj.surfaceNumVerts == 4) = 9;
    end


    function msh = getQuad4mesh(obj)
        assert(obj.cartGrid,'This method is valid only for Cartesian grids');
        msh = Mesh();
        msh.surfaces = obj.surfaces(:,1:4);
        msh.nSurfaces = size(msh.surfaces,1);
        msh.nNodes = max(msh.surfaces,[],"all");
        msh.coordinates = obj.coordinates(1:msh.nNodes,:);
        msh.surfaceNumVerts = obj.surfaceNumVerts;
        msh.surfaceNumVerts(:) = 4;
        msh.surfaceVTKType = obj.surfaceVTKType;
        msh.nDim = obj.nDim;
        %
    end

  end


  methods (Static)

    function out = copyField(fld,id)
      if ~isempty(fld)
        out = fld(id,:);
      else
        out = [];
      end

    end


    function mesh = create(varargin)

      % CREATE  Generate a mesh from an XML input struct or name-value pairs.
      %
      %   MESH = CREATE(INPUT) dispatches to the appropriate mesh constructor
      %   based on which field is present in INPUT:
      %     - 'meshFile'            : import mesh from file
      %     - 'StructuredMesh'      : build a uniform Cartesian hexahedral mesh
      %     - 'BlockStructuredMesh' : build a block-structured hexahedral mesh
      %
      % -------------------------------------------------------------------------
      % INPUTS (name-value or struct fields)
      %   meshFile             - path to an existing mesh file  [string]
      %
      %   StructuredMesh       - sub-struct with fields:
      %                            .NX, .NY, .NZ  number of cells along x, y, z
      %                            .LX, .LY, .LZ  extents [min, max]  (1x2 each)
      %
      %   BlockStructuredMesh  - sub-struct with fields:
      %                            .NX, .NY, .NZ  number of cells per block
      %                            .LX, .LY, .LZ  extents [min, max]  (1x2 each)
      %                            .depth         refinement depth (default: 2)
      %
      % -------------------------------------------------------------------------
      % OUTPUT
      %   mesh  - Mesh object built by the selected constructor
      default = struct('meshFile',missing,...
        'StructuredMesh',missing,...
        'BlockStructuredMesh',missing);

      input = readInput(default,varargin{:});


      if ~ismissing(input.meshFile)
        mesh = Mesh();
        mesh.importMesh(input.meshFile);
        return
      end

      if ~ismissing(input.StructuredMesh)
        d = double.empty;
        default = struct('NX',d,'NY',d,'NZ',d,'LX',d,'LY',d,'LZ',d);
        in = readInput(default,input.StructuredMesh);
        assert(size([in.LX;in.LY;in.LZ],2)==2,"Length of structured mesh must be a 1x2 row array")
        mesh = structuredMesh(in.NX,in.NY,in.NZ,in.LX,in.LY,in.LZ);
      elseif isfield(input,"BlockStructuredMesh")
        default = struct('NX',[],'NY',[],'NZ',[],'LX',[],'LY',[],'LZ',[],'depth',2);
        in = readInput(default,input.BlockStructuredMesh);
        assert(size([in.LX;in.LY;in.LZ],2)==2,"Length of structured mesh must be a 1x2 row array")
        mesh = BlockStructuredMesh(NX,NY,NZ,LX,LY,LZ,depth);
      end


    end
  end



end
