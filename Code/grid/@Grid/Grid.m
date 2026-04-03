classdef Grid < handle

  properties (SetAccess = public, GetAccess = public)
    
    nDim = 0    % mesh dimensions: 2(D) or 3(D)

    Cells
    Surfaces    % external faces with tag
    Faces 
    %edges      % not yet supported

    nNodes = 0
    coordinates
   
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

    % Available VTK types
    cellVTK = [10, 12, 29];     % tetra, hexa, hexa27
    surfaceVTK = [5, 9, 28];    % triangle, quad, quad9
    % edgeVTK = []

    isMixed = false             % flag for presence of multiple cell shapes
  end

  methods (Access = public)

    importMesh(obj,fileName)

    processGeometry(obj)

  
    function [surfMesh,varargout] = getSurfaceMesh(obj, surfTag, varargin)
        % Function to build a separate 2D mesh object based on the
        % surfaceTag of another object
        % there are 3 ways to call this method
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
        surfMesh = Grid();
        surfTopol = obj.Surfaces(id,:);
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

        surfMesh.surfaceNumVerts = Grid.copyField(obj.surfaceNumVerts,id);
        nNmax = max(surfMesh.surfaceNumVerts);
        surfMesh.surfaces = (reshape(surfTopol, [], nNmax));
        surfMesh.coordinates = obj.coordinates(surfOrd,:);
        surfMesh.nNodes = length(surfMesh.coordinates);
        surfMesh.nSurfaces = size(surfMesh.surfaces,1);
        surfMesh.surfaceTag = obj.surfaceTag(id);
        surfMesh.nSurfaceTag = 1;
        surfMesh.surfaceVTKType = Grid.copyField(obj.surfaceVTKType,id);
        surfMesh.nDim = 3;
        surfMesh.surfaceCentroid = Grid.copyField(obj.surfaceCentroid,id);
        surfMesh.surfaceArea = Grid.copyField(obj.surfaceArea,id);
        surfMesh.meshType = obj.meshType;
    end


    function addSurface(obj,id,topol)
       % add a surface to mesh object given the surface topology
       surf = load(topol); % standard topology file
       obj.Surfaces = [obj.Surfaces; surf];
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
        msh = Grid();
        msh.surfaces = obj.Surfaces(:,1:4);
        msh.nSurfaces = size(msh.surfaces,1);
        msh.nNodes = max(msh.surfaces,[],"all");
        msh.coordinates = obj.coordinates(1:msh.nNodes,:);
        msh.surfaceNumVerts = obj.surfaceNumVerts;
        msh.surfaceNumVerts(:) = 4;
        msh.surfaceVTKType = obj.surfaceVTKType;
        msh.nDim = obj.nDim;
        %
    end


    function nodes = getCellNodes(obj,id)
      % return a matrix of size nCells x nNode
      % expects that the input id refer to cells of the same type

      if obj.isMixed
        nodes = obj.cells.connectivity.getArray(id);
      else
        nodes = obj.cells.connectivity(id,:);
      end

    end

    function nodes = getSurfNodes(obj,id)
      % return a matrix of size nCells x nNode
      % expect that the input id refer to surfaces of the same type

      if obj.isMixed
        nodes = obj.surfaces.connectivity.getArray(id);
      else
        nodes = obj.surface.connectivity(id,:);
      end
    end


  end





  methods (Static)

    function mat = makeConnectivity(mat)

      % save space
      mat = int32(mat);

      % Possible development: support array of arrays for connectivity
      if ~all(mat(:)) && size(mat,1) > 1e6
        mat = ArrayOfArrays(mat);
      end
    end

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
        mesh = Grid();
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
