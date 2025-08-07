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

    % Function to read mesh data and store this data inside the properties
    % of the object created with the class "Mesh"
    function importGMSHmesh(obj, fileName)
      % Reading input file
      [obj.coordinates, elems, regions] = mxImportGMSHmesh(fileName);
      elems = double(elems);

      % Available types in gmsh
      cellGMSH = [4,5,12];
      surfaceGMSH = [2,3,10];

      % STORING DATA INSIDE OBJECT'S PROPERTIES
      % 3D ELEMENT DATA
      if any(obj.coordinates(:,3) ~= 0)
        obj.nDim = 3;
      else
        obj.nDim = 2;
      end
      obj.nNodes = size(obj.coordinates,1);
      % cellsID = 3D element tag for readGMSHmesh.cpp
      ID = ismember(elems(:,1), cellGMSH);
      obj.cellNumVerts = elems(ID,3);
      nVerts = max(obj.cellNumVerts);
      obj.cells = elems(ID,4:nVerts+3);
      obj.cellVTKType = mapGMSHtoVTK(obj,elems(ID,1));
      obj.cellTag = elems(ID,2);
      obj.nCells = length(obj.cellTag);
      %
      % Check for unsupported elements
      if any(~ismember(obj.cellVTKType,[10,12]))
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
      ID = ismember(elems(:,1), surfaceGMSH);
      obj.surfaceNumVerts = elems(ID,3);
      nVerts = max(obj.surfaceNumVerts);
      obj.surfaces = elems(ID,4:nVerts+3);
      obj.surfaceVTKType = mapGMSHtoVTK(obj,elems(ID,1));
      obj.surfaceTag = elems(ID,2);
      obj.nSurfaces = length(obj.surfaceTag);
      %
      % Check for unsupported elements
      if any(~ismember(obj.surfaceVTKType,[5,9]))
        error(['There are unsupported surfaces in the mesh.\n', ...
          'Supported surfaces are: - 3-node triangles       (VTKType = %d)\n', ...
          '                        - 4-node quadrilaterals  (VTKType = %d)'],5,9);
      end

      %       % 1D ELEMENT DATA
      %       % cellsID = 2D surface tag for readGMSHmesh.cpp
      %       cellsID = 1;
      %       ID = ismember(elems(:,1), cellsID);
      %       obj.edges = elems(ID,4:5);
      %       obj.edgeTag = elems(ID,2);
      %       obj.nEdges = length(obj.edgeTag);
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
      % Set all cells to tag 1 if no physical volume Tag is specified in gmsh
      if all(obj.cellTag==0)
        obj.cellTag = obj.cellTag + 1;
      end
      obj.nCellTag = max(obj.cellTag);
      obj.nSurfaceTag = max(obj.surfaceTag);
    end


    function importMesh(obj, fileName)
      % Import grid topology;
      % Admitted grid formats: vtk (unstructured grids), gmsh;
      % calling MEX functions depending on file extension
      str = split(fileName,'.');
      extension = str{2};
      switch extension
        case 'vtk'
          [obj.coordinates, elems] = mxImportVTKmesh(char(fileName));
        case 'msh'
          importGMSHmesh(obj,fileName);
          return;
      end

      elems = double(elems);

      % STORING DATA INSIDE OBJECT'S PROPERTIES
      % 3D ELEMENT DATA
      if any(obj.coordinates(:,3) ~= 0)
        obj.nDim = 3;
      else
        obj.nDim = 2;
      end
      obj.nNodes = size(obj.coordinates,1);

      ID = ismember(elems(:,1), obj.cellVTK);
      obj.cellNumVerts = elems(ID,3);
      nVerts = max(obj.cellNumVerts);
      obj.cells = elems(ID,4:nVerts+3);
      obj.cellVTKType = elems(ID,1);
      obj.cellTag = elems(ID,2);
      obj.nCells = length(obj.cellTag);
      if all(obj.cellTag==0)
        obj.cellTag = obj.cellTag + 1;
        obj.nCellTag = 1;
      else
        obj.nCellTag = max(obj.cellTag);
      end
      %
      %

      % REGIONS DATA FOR 3D ELEMENT
      if exist('regions','var')
        nRegions = length(regions);
        dims = zeros(nRegions,1);
        for i = 1 : nRegions
          dims(i) = regions(i).dim;
        end
        ID = find(dims == 3);
        for i = 1 : length(ID)
          obj.cellRegions = setfield(obj.cellRegions, regions(ID(i)).name, regions(ID(i)).ID);
        end
      end

      % 2D ELEMENT DATA
      % cellsID = 2D surface tag for readGMSHmesh.cpp

      ID = ismember(elems(:,1), obj.surfaceVTK);
      obj.surfaceNumVerts = elems(ID,3);
      nVerts = max(obj.surfaceNumVerts);
      obj.surfaces = elems(ID,4:nVerts+3);
      obj.surfaceVTKType = elems(ID,1);
      obj.surfaceTag = elems(ID,2);
      obj.nSurfaces = length(obj.surfaceTag);
      obj.nSurfaceTag = max(obj.surfaceTag);
      %
      %       %
      % REGIONS DATA FOR 2D ELEMENT
      if exist("regions",'var')
        dims = zeros(nRegions,1);
        for i = 1 : nRegions
          dims(i) = regions(i).dim;
        end
        ID = find(dims == 2);
        for i = 1 : length(ID)
          obj.surfaceRegions = setfield(obj.surfaceRegions, regions(ID(i)).name, regions(ID(i)).ID);
        end
        %
      end

      % 1D ELEMENT DATA
      % cellsID = 2D surface tag for readGMSHmesh.cpp
      ID = ismember(elems(:,1), obj.edgeVTK);
      obj.edgeNumVerts = elems(ID,3);
      nVerts = max(obj.edgeNumVerts);
      obj.edges = elems(ID,4:nVerts+3);
      obj.edgeTag = elems(ID,2);
      obj.edgeVTKType = elems(ID,1);
      obj.nEdges = length(obj.edgeTag);
    end


    
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


    function createCartesianGrid(obj,dim,deg,varargin)
        % Generate Cartesian mesh of quadrilateral elements
        assert(isempty(obj.coordinates),['A mesh hase been' ...
            'already defined for this object istance'])
        switch dim
            case 2
                assert(nargin==7,['Incorrect number of input arguments for' ...
                    '2D cartesian mesh'])
                [x,y,nx,ny] = deal(varargin{1},varargin{2},varargin{3},varargin{4});
                genCartGrid2D(obj,deg,x,y,nx,ny);
            case 3
                assert(nargin==9,['Incorrect number of input arguments for ' ...
                    '3D cartesian mesh'])
                genCartGrid3D(obj,deg,varargin);
        end
        obj.cartGrid = true;
    end

    function setCartGridFace(obj,faceTag,id)      
        assert(obj.cartGrid, 'The mesh object is not a Cartesian Grid');
        % Assign a tag to specific edges/surfaces for a cartesian grid
        assert(length(faceTag) == length(id), "Number of faces does not" + ...
            "match number of id" );
        if obj.nDim < 3
            fStr = ['north','south','west','east'];
            assert(any(strcmpi(fStr,faceTag)),['Undefined face Tag: ' ...
                'Admitted tags are: north,south,west,east']);

        else
            fStr = ['top','bottom','north','south','west','east'];
                        assert(any(strcmpi(fStr,faceTag)),['Undefined face Tag: ' ...
                            'Admitted tags are: top,bottom,north,south,west,east']);
        end
        [x,y,z] = deal(obj.coordinates(:,1),obj.coordinates(:,2),obj.coordinates(:,3));
        % Retrieve list of nodes of the face based on the coordinates 
        for i = 1:length(id)
            switch faceTag(i)
                case 'north'
                    nodeList = sort(find(y == max(y)));
                case 'south'
                    nodeList = sort(find(y == min(y)));
                case 'east'
                    nodeList = sort(find(x == max(x)));
                case 'west'
                    nodeList = sort(find(x == min(x)));
                case 'top'
                    nodeList = sort(find(z == max(z)));
                case 'bottom'
                    nodeList = sort(find(z == min(z)));
            end

            % build face topology
            if obj.nDim < 3
                topol = [nodeList(1) repelem(nodeList(2:end-1),2), nodeList(end)];
                topol = (reshape(topol, 2, []))';
                obj.edges = [obj.edges; topol];
                obj.edgeTag = [obj.edgeTag id(i)*ones(size(obj.edges,1))];
            else
                 
            end

        end

    end

    

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
    end

%     function centroids = getCellCentroids(obj)
%        centroids = zeros(obj.nCells,3);
%        for i = 1:obj.nCells
%           coord = obj.coordinates(obj.cells(i,:),:);
%           centroids(i,:) = sum(coord,1)/size(coord,1);
%        end
%     end
% 
%     function centroids = getSurfaceCentroids(obj)
%       centroids = zeros(obj.nSurfaces,3);
%       for i = 1:obj.nSurfaces
%         coord = obj.coordinates(obj.surfaces(i,:),:);
%         centroids(i,:) = sum(coord,1)/size(coord,1);
%       end
%     end

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

  methods (Access = private)
      function obj = genCartGrid2D(obj,deg,x,y,nx,ny)
          % Node coordinates
          switch deg
              case 1
                  nNodElem = 4;
                  xc = linspace(x(1),x(2),nx+1);
                  yc = linspace(y(1),y(2),ny+1);
                  [X,Y] = meshgrid(xc,yc);
                  X = X'; Y=Y';
                  obj.nNodes = length(xc)*length(yc);
                  obj.coordinates = zeros(obj.nNodes,3);
                  obj.coordinates(:,1:2) = [X(:) Y(:)];
                  % Topology
                  s1 = [1 2 length(xc)+2 length(xc)+1];
                  s2 = reshape(0:((nx+1)*ny-1),nx+1,[]);
                  s2 = s2(1:end-1,:);
                  obj.surfaces = s1+s2(:);
              case 2
                  nNodElem = 8;
                  xc = linspace(x(1),x(2),nx+1);
                  yc = linspace(y(1),y(2),ny+1);
                  [X,Y] = meshgrid(xc,yc);
                  X = X'; Y=Y';
                  coord1 = [X(:) Y(:)]; % grid of angle nodes
                  dx = 0.5*(x(2)-x(1))/nx;
                  dy = 0.5*(y(2)-y(1))/ny;
                  xc = linspace(x(1)+dx,x(2)-dx,nx);
                  yc = linspace(y(1),y(2),ny+1);
                  [X,Y] = meshgrid(xc,yc); X = X'; Y=Y';
                  coord2 = [X(:) Y(:)]; % grid of horizontal edge nodes
                  xc = linspace(x(1),x(2),nx+1);
                  yc = linspace(y(1)+dy,y(2)-dy,ny);
                  [X,Y] = meshgrid(xc,yc); X = X'; Y=Y';
                  coord3 = [X(:) Y(:)]; % grid of vertical edge nodes
                  obj.nNodes = size(coord1,1)+size(coord2,1)+size(coord3,1);
                  obj.coordinates = zeros(obj.nNodes,3);
                  obj.coordinates(:,1:2) = [coord1;coord2;coord3];
                  % TOPOLOGY
                  k = 0;
                  top = zeros(nx*ny,8);
                  for iy = 1:ny
                      for ix = 1:nx
                          k = k+1;
                          n = ix+(iy-1)*(nx+1); 
                          n1 = [n n+1 n+nx+2 n+nx+1]; %nodes 1 2 3 4
                          n2 = size(coord1,1)+[ix+(iy-1)*nx ix+iy*nx]; %nodes 5 7
                          n3 = size(coord1,1)+size(coord2,1)+...
                              [ix+(iy-1)*(nx+1)+1 ix+(iy-1)*(nx+1)]; %nodes 6 8
                          top(k,:) = [n1 n2(1) n3(1) n2(2) n3(2)];
                      end
                  end
                  obj.surfaces = top;
          end

          if any(obj.coordinates(:,3) ~= 0)
              obj.nDim = 3;
          else
              obj.nDim = 2;
          end

          obj.nSurfaces = size(obj.surfaces,1);
          obj.surfaceNumVerts = nNodElem*ones(obj.nSurfaces,1);
          obj.surfaceVTKType = 9*ones(obj.nSurfaces,1);
          obj.surfaceTag = ones(obj.nSurfaces,1); 
          % surfaceTag property can be modified using specifc method for
          % CartGrids

          % 1D ELEMENT DATA
          % cellsID = 2D surface tag for readGMSHmesh.cpp
          obj.nEdges = 0;
          % Edge datas for CartGrid are introduced using the setGridFace method 
          obj.nSurfaceTag = max(obj.surfaceTag);
      end

      function outList = mapGMSHtoVTK(obj,list)
        % map a list of gmsh id to corresponding vtk id
        gmshId = [2,3,10,4,5,12];
        vtkId = [obj.surfaceVTK obj.cellVTK];
        outList = zeros(length(list),1);
        for i = 1:length(gmshId)
          id = list == gmshId(i);
          outList(id) = vtkId(i);
        end
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
  end



end
