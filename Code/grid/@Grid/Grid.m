classdef Grid < handle

  properties (SetAccess = public, GetAccess = public)
    
    nDim = 0    % mesh dimensions: 2(D) or 3(D)

    cells = struct('connectivity',[]); 
    surfaces = struct('connectivity',[]);    % external faces with tag
    faces 
    edges         

    nNodes = 0
    coordinates
   
  end

  properties (GetAccess = private, SetAccess = private)


    % edgeVTK = []

    isProcessed = false

    isMixed = false             % flag for presence of multiple cell shapes
    structFlag = false        % flag for structured mesh
  end

  methods (Access = public)

    importMesh(obj,fileName)

    initializeGrid(obj)

    processGeometry(obj)

    [vol,center] = computeCellGeometry(obj)

    processFaces(obj,vtkId)

    processEdges(obj,vtkId)

    function obj = Grid()
      initializeGrid(obj)
    end


    function outGrid = getSurfaceGrid(obj, surfTag)

      surf = obj.surfaces;
      outGrid = Grid();

      % select surfaces
      if nargin == 1
        id = true(surf.num,1);
      elseif nargin == 2
        if islogical(surfTag)
          assert(numel(surfTag) == surf.num);
          id = surfTag(:);
        else
          id = ismember(surf.tag, surfTag);
        end
      else
        error("getSurfaceGrid:wrong number of input arguments")
      end


      if ~any(id); return; end

      sel  = find(id);
      nSel = numel(sel);

      % extract connectivity
      connSel = getRows(surf.connectivity, sel);

      if isa(connSel,'ArrayOfArrays')
        [oldConn, ~] = connSel.getData();
        nVert = connSel.arraySize();
        nVert = nVert(:);
      else
        oldConn = connSel.';
        oldConn = oldConn(:);
        nVert = repmat(size(connSel,2), nSel, 1);
      end

      % stable local node renumbering
      [globalNodeId, ~, newConn] = unique(oldConn, 'stable');

      % local connectivity
      if all(nVert == nVert(1))
        surfConn = reshape(newConn, nVert(1), []).';
      else
        surfConn = ArrayOfArrays(newConn, nVert);
      end

      % populate output grid
      surfOut = struct();
      surfOut.connectivity = surfConn;
      surfOut.numVerts     = nVert;
      surfOut.num          = nSel;
      surfOut.tag          = surf.tag(sel);
      surfOut.nTag         = numel(unique(surfOut.tag));
      surfOut.VTKType      = surf.VTKType(sel);
      surfOut.vtkTypes     = reshape(unique(surfOut.VTKType),1,[]);
      surfOut.center       = surf.center(sel,:);
      surfOut.normal       = surf.normal(sel,:);
      surfOut.area         = surf.area(sel);
      surfOut.faceId       = surf.faceId(sel);

      edgeStruct.num = 0;
      edgeStruct.length = zeros(0,1);
      edgeStruct.center = zeros(0,3);
      edgeStruct.neighbors = zeros(0,2);    % edge to surface
      edgeStruct.isBoundary = false(0,1);
      edgeStruct.connectivity = zeros(0,2); % edge to node

      % local 2 global node map
      surfOut.loc2glob = globalNodeId(:);
      surfOut.surfaces2edges = ArrayOfArrays();

      outGrid.coordinates = obj.coordinates(globalNodeId,:);
      outGrid.nNodes      = size(outGrid.coordinates,1);
      outGrid.nDim        = 3;
      outGrid.surfaces    = surfOut;
      outGrid.edges = edgeStruct;

    end
    %
    %
    % function addSurface(obj,id,topol)
    %    % add a surface to mesh object given the surface topology
    %    surf = load(topol); % standard topology file
    %    obj.surfaces = [obj.surfaces; surf];
    %    if any(obj.surfaceTag==id)
    %       error('Surface iD already been defined')
    %    else
    %       obj.surfaceTag = [obj.surfaceTag; id*ones(size(surf,1),1)];
    %    end
    % 
    %    if isempty(obj.nSurfaceTag)
    %       obj.nSurfaceTag = 1;
    %    else
    %       obj.nSurfaceTag = obj.nSurfaceTag + 1;
    %    end
    %    obj.surfaceNumVerts = [obj.surfaceNumVerts; sum(surf > 0,2)];
    %    obj.surfaceVTKType(obj.surfaceNumVerts == 3) = 5;
    %    obj.surfaceVTKType(obj.surfaceNumVerts == 4) = 9;
    % end
    % 
    % 
    % function msh = getQuad4mesh(obj)
    %     assert(obj.cartGrid,'This method is valid only for Cartesian grids');
    %     msh = Grid();
    %     msh.surfaces = obj.surfaces(:,1:4);
    %     msh.nSurfaces = size(msh.surfaces,1);
    %     msh.nNodes = max(msh.surfaces,[],"all");
    %     msh.coordinates = obj.coordinates(1:msh.nNodes,:);
    %     msh.surfaceNumVerts = obj.surfaceNumVerts;
    %     msh.surfaceNumVerts(:) = 4;
    %     msh.surfaceVTKType = obj.surfaceVTKType;
    %     msh.nDim = obj.nDim;
    %     %
    % end




    function outList = getCellsByVTKId(obj,vtkId,cellList)
      if ~isscalar(vtkId)
        error("Input vtk id must be a scalar value");
      end

      if nargin == 2
        cellList = (1:obj.cells.num)';
      end

      outList = cellList(obj.cells.VTKType(cellList) == vtkId);
    end

    function surfList = getSurfByVTKId(obj,vtkId)
      if ~isscalar(vtkId)
        error("Input vtk id must be a scalar value");
      end
      surfList = find(obj.surfaces.VTKType == vtkId);
    end


    function nodes = getCellNodes(obj,id)
      % return a matrix of size nCells x nNode
      % expects that the input id refer to cells of the same type

      if nargin == 1
        id = 1:obj.cells.num;
      end

      nodes = getRowsMatrix(obj.cells.connectivity,id);

    end

    function nodes = getSurfNodes(obj,id)
      % return a matrix of size nCells x nNode
      % expect that the input id refer to surfaces of the same type


      if nargin == 1
        id = 1:obj.surfaces.num;
      end

       nodes = getRowsMatrix(obj.surfaces.connectivity,id);

    end


    function setConnectivity(obj,type,rows,conn)
      % modify a chunk of the connectivity matrix of elements of the same size

      if ~any(strcmp(type,["cells","surfaces"]))
        error("Invalid connectivity type specifier. Acceptet values are 'cells' or 'surfaces'")
      end

      if obj.isMixed
        obj.(type).connectivity.setArray(rows,conn);
      else
        obj.(type).connectivity(rows,:) = conn;
      end
    end


    function setStructured(obj)
      obj.structFlag = true;
    end

    function out = isGeometryProcessed(obj)

      out = obj.isProcessed;
      
    end


    function out = isStructured(obj)
      out = obj.structFlag;
    end

    function computeAvgNodalNormal(obj)

      avNorm = zeros(obj.nNodes,3);
      coord = obj.coordinates;
      surf = obj.surfaces;
      topol = surf.connectivity;

      for i = 1:surf.num
        nVert = surf.numVerts(i);
        coordLoc = coord(topol(i,[nVert,1:nVert,1]),:);
        % get edges direction according to predefined node ordering
        dn = diff(coordLoc,1);
        for j = 1:nVert
          nLoc = cross(dn(j,:),dn(j+1,:));
          avNorm(topol(i,j),:) = avNorm(topol(i,j),:) + nLoc;
        end
      end
      % normalize matrix
      normN = sqrt(sum(avNorm.^2,2));
      obj.surfaces.avgNodNormal = avNorm./(normN);
    end



    function [list,infl] = getNodeInfluence(obj,srcEnt,entId,fl)

      % return the nodal incidence of a subset of source entities
      % list is an array of array with numel(entId) subarrays
      % infl is the influence value for each node in each cell
      % if fl == "interp", the influence is interpolative 

      % process by vtk type
      switch srcEnt
        case entityField.cell
          mesh = obj.cells;
          getNodes = @(el) obj.getCellNodes(el);
        case entityField.surface
          mesh = obj.surfaces;
          getNodes = @(el) obj.getSurfNodes(el);
      end

      interp = false;
      if nargin > 3 && strcmpi(fl,"interp")
        interp = true;
      end

      vtkTypes = mesh.vtkTypes;
      coords = obj.coordinates;

      list = cell(length(vtkTypes),1);
      infl = cell(length(vtkTypes),1);

      for i = 1:length(vtkTypes)
        vtkId = vtkTypes(i);

        elem = FiniteElementType.create(vtkId,obj);
        elIds = find(mesh.VTKType == vtkId);

        elIds = intersect(entId,elIds);

        topol = getNodes(elIds);

        listLoc = getIncidenceID(entityField.node,obj,srcEnt,elIds);
        infLoc = zeros(listLoc.totalSize,1);

        list{i} = listLoc;

        % element type for processing
        k = 0;

        for j = 1:length(elIds)
          % nodes = topol(j,:);
          % coord = coords(nodes,:);
          infEl = elem.getNodeInfluence(elIds(j));
          if interp
            infEl = infEl/sum(infEl);
          end
          nn = listLoc.arraySize(j);
          infLoc(k+1:k+nn) = infEl;

          k = k+nn;
        end

        infl{i} = infLoc;

      end

      list = vertcat(list{:});
      infl = cell2mat(infl);

    end


    function list = getFlatConnectivity(obj,fld)

      % str: either cells or surfaces
      conn = obj.(fld).connectivity;

      if isa(conn,"ArrayOfArrays")
        % connectivity is an ArrayOfArrays
        list = conn.getData;
      else % connectivity is a strandard matrix
        conn = conn';
        list = conn(:);
      end

    end

  end





  methods (Static)


    function mat = makeConnectivity(mat)

      % save space
      % mat = int32(mat);

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


    function grid = create(varargin)

      % CREATE  Generate a mesh from an XML input struct or name-value pairs.
      %
      %   GRID = CREATE(INPUT) dispatches to the appropriate mesh constructor
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
        grid = Grid();
        grid.importMesh(input.meshFile);
        return
      end

      if ~ismissing(input.StructuredMesh)
        d = double.empty;
        default = struct('NX',d,'NY',d,'NZ',d,'LX',d,'LY',d,'LZ',d);
        in = readInput(default,input.StructuredMesh);
        assert(size([in.LX;in.LY;in.LZ],2)==2,"Length of structured mesh must be a 1x2 row array")
        grid = structuredMesh(in.NX,in.NY,in.NZ,in.LX,in.LY,in.LZ);
      elseif isfield(input,"BlockStructuredMesh")
        default = struct('NX',[],'NY',[],'NZ',[],'LX',[],'LY',[],'LZ',[],'depth',2);
        in = readInput(default,input.BlockStructuredMesh);
        assert(size([in.LX;in.LY;in.LZ],2)==2,"Length of structured mesh must be a 1x2 row array")
        grid = BlockStructuredMesh(NX,NY,NZ,LX,LY,LZ,depth);
      end


    end
  end



end
