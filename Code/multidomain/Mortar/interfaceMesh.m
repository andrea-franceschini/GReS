classdef interfaceMesh < handle
  % Manager for topological information and operation on interfaces between
  % meshes

  properties (Access = public)
    activeCells
  end

  
  properties
    % prop{1} -> master     prop{2} -> slave 
    msh
    local2glob
    elemConnectivity
    nN
    nEl
    cellType
    nEdges
    e2n
    n2e
    e2f
    f2c
    f2e
    avgNodNormal    % average nodal normal for master and slave side
  end
  
  methods
    function obj = interfaceMesh(domains,surf)
      mshMaster = domains(1).grid.topology;
      mshSlave = domains(2).grid.topology;
      % extrad 2d surface mesh from parent domain
      obj.msh = [mshMaster.getSurfaceMesh(surf{1});
                 mshSlave.getSurfaceMesh(surf{2})];
      getCellTypes(obj);
      getConnectivityMatrix(obj);
      mapLocalNod2Glob(obj,1,mshMaster,surf{1});
      mapLocalNod2Glob(obj,2,mshSlave,surf{2});
      setupEdgeTopology(obj,1);
      setupEdgeTopology(obj,2);
      obj.buildFace2CellMap([mshMaster mshSlave]);
      obj.computeAverageNodalNormal();
    end

    function getConnectivityMatrix(obj)
      cs = ContactSearching(obj.msh(1),obj.msh(2));
      obj.elemConnectivity = cs.elemConnectivity;
      obj.activeCells{2} = find(any(obj.elemConnectivity,1));
      obj.activeCells{1} = reshape(find(any(obj.elemConnectivity,2)),1,[]);
      obj.nEl(2) = numel(obj.activeCells{2});
      obj.nEl(1) = sum(any(obj.elemConnectivity,2));
    end

    function list = getActiveCells(obj,side,id)
      % side: 1 -> master side:2 -> slave
      if nargin == 2
        list = obj.activeCells{side};
      elseif nargin == 3
        list = obj.activeCells{side}(id);
      else
        error('Too many input arguments')
      end
    end
  end


  methods (Access = private)

    function setupEdgeTopology(obj, sideID)
      % Extended to handle both triangular and quadrilateral faces
      % inspired by face topology for FV in MRST

      mesh = obj.msh(sideID);
      faces = mesh.surfaces;
      nF = mesh.nSurfaces;


      switch obj.cellType(sideID)
        case 10
          % Triangular elements: edges are [1,2], [2,3], [3,1]
          e1 = faces(:, [1, 2]);  % edge 1
          e2 = faces(:, [2, 3]);  % edge 2
          e3 = faces(:, [3, 1]);  % edge 3

          edgeMat = [e1; e2; e3];

          % Edge local indices
          id = [ (1:nF)', repmat(1, nF, 1);    % edge 1
            (1:nF)', repmat(2, nF, 1);    % edge 2
            (1:nF)', repmat(3, nF, 1) ];  % edge 3

          nEdgesPerFace = 3;

        case 12
          % Quadrilateral elements: edges are [1,2], [2,3], [3,4], [4,1]
          bot = faces(:, [1, 2]);  % bottom edge
          rgt = faces(:, [2, 3]);  % right edge
          top = faces(:, [3, 4]);  % top edge
          lft = faces(:, [4, 1]);  % left edge

          edgeMat = [bot; rgt; top; lft];

          id = [ (1:nF)', repmat(1, nF, 1);    % bottom
            (1:nF)', repmat(2, nF, 1);    % right
            (1:nF)', repmat(3, nF, 1);    % top
            (1:nF)', repmat(4, nF, 1) ];  % left

          nEdgesPerFace = 4;

        otherwise
          error('Unsupported element type: each face must have 3 or 4 vertices.');
      end

      % Sort nodes in each edge (edge = {min,max})
      [edgeMat, i] = sort(edgeMat, 2);
      i = i(:,1);

      % Sort edgeMat rows to group duplicates
      [edgeMat, j] = sortrows(edgeMat);

      % Run-length encode to get unique edges and their multiplicity
      [obj.e2n{sideID}, n] = rle(edgeMat);
      obj.nEdges(sideID) = numel(n);

      % Map each occurrence back to edge index
      N = repelem(1:obj.nEdges(sideID), n);  % same size as edgeMat
      edgeIDs = zeros(size(edgeMat, 1), 1);
      edgeIDs(j) = N;

      % Reshape to f2e mapping: nF rows, 3 or 4 edges per face
      obj.f2e{sideID} = reshape(edgeIDs, [nF, nEdgesPerFace]);

      % Build e2f mapping (each edge can be shared by up to 2 faces)
      obj.e2f{sideID} = accumarray([N', i(j)], id(j,1), [obj.nEdges(sideID), 2]);

      allEdges = obj.e2n{sideID};
      maxNode = max(allEdges(:));
      obj.n2e{sideID} = cell(maxNode, 1);
      for e = 1:size(allEdges,1)
        nodes = allEdges(e,:);
        obj.n2e{sideID}{nodes(1)}(end+1) = e;
        obj.n2e{sideID}{nodes(2)}(end+1) = e;
      end
    end

    % provisional: this method will be removed to support different element
    % types in each side. 
    function getCellTypes(obj)
      for i = [1 2]
        switch obj.msh(i).surfaceVTKType(1)
          case 5 % Triangle mesh Master
            obj.cellType(i) = 10;
            obj.nN(i) = 3;
          case {9,28} % Quad mesh Master
            obj.cellType(i) = 12;
            obj.nN(i) = 4;
        end
      end
    end


    function mapLocalNod2Glob(obj,side,msh,surf)
      % return array where local node numbering map to node id in the full
      % 3D grid
      surfGlob2loc = find(ismember(msh.surfaceTag,surf));
      globNodes = (msh.surfaces(surfGlob2loc,:))';
      globNodes = globNodes(:);
      obj.local2glob{side} = zeros(obj.msh(side).nNodes,1);
      locNodes = (obj.msh(side).surfaces)';
      locNodes = locNodes(:);
      obj.local2glob{side}(locNodes) = globNodes;
    end

    function buildFace2CellMap(obj, meshBg)

      for i = [1 2]
        top = obj.local2glob{i}(obj.msh(i).surfaces);  % global face node IDs
        nFaces = obj.nEl(i);
        nFaceNodes = size(top, 2);

        % Initialize face-to-cell mapping
        obj.f2c{i} = zeros(nFaces, 1);

        % Build node-to-cell adjacency list
        allCells = meshBg(i).cells;  % size: [nCells, nodesPerCell]
        nCells = size(allCells, 1);
        maxNodeID = max(allCells(:));

        node2cells = cell(maxNodeID, 1);  % preallocate
        for c = 1:nCells
          for v = allCells(c, :)
            node2cells{v} = [node2cells{v}, c];
          end
        end

        % For each face, find the unique cell that contains all its nodes
        for e = 1:nFaces
          faceNodes = top(e, :);

          % Get candidate cells as the intersection of node2cell lists
          candidates = node2cells{faceNodes(1)};
          for k = 2:nFaceNodes
            candidates = intersect(candidates, node2cells{faceNodes(k)});
            if isempty(candidates)
              break;
            end
          end

          % Among candidates, find the one that contains all the face nodes
          hasCellNeigh = false;
          for id = candidates
            if all(ismember(faceNodes, allCells(id, :)))
              obj.f2c{i}(e) = id;
              hasCellNeigh = true;
              break;
            end
          end

          assert(hasCellNeigh, 'Invalid connectivity: face does not belong to any cell.');
        end
      end
    end

    function computeAverageNodalNormal(obj)
      for side = 1:2
        avNorm = zeros(obj.msh(side).nNodes,3);
        coord = obj.msh(side).coordinates;
        topol = obj.msh(side).surfaces;
        nv = obj.msh(side).surfaceNumVerts;
        for el = 1:obj.msh(side).nSurfaces
          nVert = nv(el);
          coordLoc = coord(topol(el,[nVert,1:nVert,1]),:);
          % get edges direction according to predefined node ordering
          dn = diff(coordLoc,1);
          for i = 1:nv(el)
            nLoc = cross(dn(i,:),dn(i+1,:));
            avNorm(topol(el,i),:) = avNorm(topol(el,i),:) + nLoc;
          end
        end
        % normalize matrix
        normN = sqrt(sum(avNorm.^2,2));
        obj.avgNodNormal{side} = avNorm./(normN);
      end
    end


  end

end

