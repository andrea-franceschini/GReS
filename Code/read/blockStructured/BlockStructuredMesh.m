classdef BlockStructuredMesh < handle
  % implements the logic of a blockstructured mesh with octree logic
  % refinement

  properties
    blocks = Block.empty % 3D array of blocks trees
    dims
    nCells
    maxDepth
    mapDiffToFace = [2,3,1,4,5,6]
    faceOpposite = [6,4,5,2,3,1]
  end

  methods
    function obj = BlockStructuredMesh(dims,Ncells,maxDepth)
      % dims: 3x2 with extreme coordinate of box grid
      % nCells: 1x3 array with number of cell for each direction
      % maxDepth: maximum depth of recursive refinement
      obj.dims = dims;
      obj.nCells = Ncells;
      obj.maxDepth = maxDepth;
      createBaseGrid(obj)
    end

    function createBaseGrid(obj)
      obj.blocks(obj.nCells(1),obj.nCells(2),obj.nCells(3)) = Block();
      xVec = linspace(obj.dims(1,1),obj.dims(1,2),obj.nCells(1)+1);
      yVec = linspace(obj.dims(2,1),obj.dims(2,2),obj.nCells(2)+1);
      zVec = linspace(obj.dims(3,1),obj.dims(3,2),obj.nCells(3)+1);
      for k = 1:obj.nCells(3)
        for j = 1:obj.nCells(2)
          for i = 1:obj.nCells(1)
            bb = [xVec(i) xVec(i+1); yVec(j) yVec(j+1); zVec(k) zVec(k+1)];
            obj.blocks(i,j,k) = Block(bb,0,[0,0,0]);
          end
        end
      end
    end

    function refineGrid(obj,targetPoint,varargin)
      % get id of block containing the target point
      xVec = linspace(obj.dims(1,1),obj.dims(1,2),obj.nCells(1)+1);
      yVec = linspace(obj.dims(2,1),obj.dims(2,2),obj.nCells(2)+1);
      zVec = linspace(obj.dims(3,1),obj.dims(3,2),obj.nCells(3)+1);
      i = sum(xVec < targetPoint(1));
      j = sum(yVec < targetPoint(2));
      k = sum(zVec < targetPoint(3));
      if isempty(varargin)
        md = obj.maxDepth;
      else
        md = varargin{1};
      end
      refineBlock(obj.blocks(i,j,k),md,targetPoint);
    end

    function refineRecursive(obj,blockPosition,maxDepth,varargin)
      % varargin: an array for lower level refinement
      % each entry tells the lower level children target
      % children are numbered from 1 to 8 according to the structured
      % indexing rule
      i = blockPosition(1);
      j = blockPosition(2);
      k = blockPosition(3);
      if isempty(varargin)
        refineBlock(obj.blocks(i,j,k),maxDepth);
      else
        child = obj.blocks(i,j,k);
        lev = varargin{1};
        for i = 1:numel(varargin{1})
          child = child.children(lev(i));
        end
        refineBlock(child,maxDepth);
      end
    end 

    function [leaves,leavesBlocks] = getLeaves(obj)

      % get leaves and store associated global indices
      maxLeaves = prod(obj.nCells)*8^(obj.maxDepth);
      leaves(maxLeaves,1) = Block();
      leavesBlocks = zeros(maxLeaves,3);
      l=0;
      for k = 1:obj.nCells(3)
        for j = 1:obj.nCells(2)
          for i = 1:obj.nCells(1)
            leave_i = getLeaves(obj.blocks(i,j,k));
            nl = numel(leave_i);
            leaves(l+1:l+nl) = leave_i;
            leavesBlocks(l+1:l+nl,:) = repmat([i,j,k],nl,1);
            l = l+nl;
          end
        end
      end

      leaves = leaves(1:l);
      leavesBlocks = leavesBlocks(1:l,:);

    end

    function mesh = processGeometry(obj)
      % build a GReS Mesh() based on block structured grid
      % loop over leaves and process adjacency
      % mark each face of the leaves as conforming, master, slave
      % for each leaf, add new independent nodes to the list
      % add topology of non conforming faces (slave or master depending on
      % resolution)
      % remove duplicated nodes on conforming faces

      [leaves, leavesBlocks] = getLeaves(obj);

      cellOrder = zeros(numel(leaves),1);
      coordinates = zeros(8*numel(leaves),3);

      nodeNeighbors = (1:size(coordinates,1))';

      k=0;

      % raw topology with 8 distinct nodes for each cell
      cells = (1:numel(nodeNeighbors))';
      surfaces = zeros(24*size(cells,1),1);
      surfaceTag = zeros(size(surfaces,1),1);

      for i = 1:numel(leaves)
        % define all nodes coordinates (DG-like)
        coordinates(8*(i-1)+1:8*(i-1)+8,:) = getNodeCoordinates(leaves(i));
      end


      ns1 = 0;
      ns2 = 0;
      for i = 1:size(leavesBlocks,1)
        diff = abs(leavesBlocks(i+1:end,:) - leavesBlocks(i,:));
        neighbors = i+find(sum(diff,2) < 2);
        for j = neighbors'

          [faceId, faceTag] = compareCells(obj,leaves(i),leaves(j),leavesBlocks(i,:),leavesBlocks(j,:));
          % mark adjacent faces as conforming, master or slave
          if isempty(faceId)
            % non neighbors
            continue
          elseif leaves(i).faceTag(faceId(1))>=0 && leaves(j).faceTag(faceId(2))>=0
            % interface already processed
            continue
          end

          id = [i,j];

          leaves(i).faceTag(faceId(1)) = faceTag(1);
          leaves(j).faceTag(faceId(2)) = faceTag(2);


          % update surfaces and topology depending on the tag
          if all(faceTag == 0)
            % CONFORMING FACE
            if cellOrder(i)==0
              k = k+1;
              cellOrder(i) = k;
            end
            if cellOrder(j)==0
              k = k+1;
              cellOrder(j) = k+1;
            end
            % check leading cell
            [~,m1] = min(cellOrder([i,j]));
            node = {8*(id(m1)-1) + Block.faceToNode(faceId(m1),:);
              8*(id(3-m1)-1) + Block.faceToNode(faceId(3-m1),:)};
            % update node neighbors from leading cell
            nodeNeighbors(node{2}) = nodeNeighbors(node{1});
          else
            % NON CONFORMING FACE
            % get global index of nodes of adjacent faces
            [~,m2] = max([leaves(i).level, leaves(j).level]);
            surfaces(ns2+1:ns2+4) = 8*(id(m2)-1) + Block.faceToNode(faceId(m2),:);
            surfaces(ns2+5:ns2+8) = 8*(id(3-m2)-1) + Block.faceToNode(faceId(3-m2),:);
            surfaceTag(ns1+1) = 1;
            surfaceTag(ns1+2) = 2;
            ns1 = ns1 + 2;
            ns2 = ns2 + 8;
          end
        end
      end

      surfaces = surfaces(1:ns2);
      surfaceTag = surfaceTag(1:ns1);

      % remove duplicated nodes in conforming faces
      [~,in,ic] = unique(nodeNeighbors);
      %nodeNeighborsUnique = nodeNeighbors(ic);
      surfaces = ic(surfaces);
      cells = ic(cells);
      %surfaces = nodeNeighborsUnique(surfaces);
      cells = (reshape(cells,8,[]))';
      surfaces = (reshape(surfaces,4,[]))';
      [~,is,~] = unique(sort(surfaces,2),'rows','stable');
      surfaces = surfaces(is,:);
      surfaceTag = surfaceTag(is);
      coordinates = coordinates(nodeNeighbors(in),:);

      % add external surfaces for fast bcs imposition
      % 3: bottom
      % 4: top
      % 5: south
      % 6: north
      % 7: east
      % 8: west

      allFaces = [cells(:,[1,2,3,4]);
                  cells(:,[5,6,7,8]);
                  cells(:,[1,2,6,5]);
                  cells(:,[4,3,7,8]);
                  cells(:,[1,4,8,5]);
                  cells(:,[2,3,7,6])];

      xMin = min(coordinates(:,1));
      xMax = max(coordinates(:,1));
      yMin = min(coordinates(:,2));
      yMax = max(coordinates(:,2));
      zMin = min(coordinates(:,3));
      zMax = max(coordinates(:,3));

      x = reshape(coordinates(allFaces,1),[],4);
      y = reshape(coordinates(allFaces,2),[],4);
      z = reshape(coordinates(allFaces,3),[],4);

      isBottom = all(z==zMin,2);
      isTop = all(z==zMax,2);
      isSouth = all(y==yMin,2);
      isNorth = all(y==yMax,2);
      isEast = all(x==xMin,2);
      isWest = all(x==xMax,2);

      N = [sum(isBottom); sum(isTop); 
           sum(isSouth); sum(isNorth);
           sum(isEast); sum(isWest)];

      surfaces = [surfaces;
                 [allFaces(isBottom,:);
                  allFaces(isTop,:);
                  allFaces(isSouth,:);
                  allFaces(isNorth,:);
                  allFaces(isEast,:);
                  allFaces(isWest,:)]];

      if ~isempty(surfaceTag)
        s = max(surfaceTag);
      else
        s = 0;
      end

      surfaceTag = [surfaceTag;
                    repelem((s+1:s+6)',N)];


              %        v
    % 4----------3
    % |\     ^   |\
    % | \    |   | \
    % |  \   |   |  \
    % |   8------+---7
    % |   |  +-- |-- | -> u
    % 1---+---\--2   |
    %  \  |    \  \  |
    %   \ |     \  \ |
    %    \|      w  \|
    %     5----------6



      


      % populate Mesh object
      mesh = Mesh();
      mesh.coordinates = coordinates;
      mesh.nNodes = size(coordinates,1);
      mesh.cells = cells;
      mesh.nCells = size(cells,1);
      mesh.cellTag = ones(mesh.nCells,1);
      mesh.nSurfaces = length(surfaceTag);
      mesh.surfaceTag = surfaceTag;
      mesh.surfaces = surfaces;
      mesh.nCellTag = max(mesh.cellTag);
      mesh.nSurfaceTag = max(mesh.surfaceTag);
      mesh.cellVTKType = 12*ones(mesh.nCells,1);
      mesh.surfaceVTKType = 9*ones(mesh.nSurfaces,1);
      mesh.cellNumVerts = 8*ones(mesh.nCells,1);
      mesh.surfaceNumVerts = 4*ones(mesh.nSurfaces,1);
      mesh.nDim = 3;

    end


    function [faceId,faceTag] = compareCells(obj,block1,block2,rootId1,rootId2)
      % check the level
      % if same, possibly conforming neighbors
      % if different, non conforming neighbors, higher level is master,
      % lower level is slave
      % to check if they are neighbors, get a common index in the finer
      % level grid. check adjacency along a direction and overlap on the
      % other two

      % get global index of each block
      levels = [block1.level, block2.level];

      % get global indices for adjacency comparison

      [~,m] = max(levels);

      % scale to finer grid
      i1 = block1.globalIndex*2^(levels(m)-levels(1));
      i2 = block2.globalIndex*2^(levels(m)-levels(2));

      % compare in global mesh numbering
      % compare global cell index in finer grid
      i1 = (rootId1-1)*2^max(levels) + i1 + [0, 2^(levels(m)-levels(1))-1]';
      i2 = (rootId2-1)*2^max(levels) + i2 + [0, 2^(levels(m)-levels(2))-1]';

      faceId1 = obj.getAdjacentFace(i1,i2);

      if isempty(faceId1)
        % same cell or non neighbor cell
        faceId = []; faceTag = [];
        return
      else
        % check neighboring face
        faceId2 = obj.faceOpposite(faceId1);
        if levels(1) == levels(2)
          % mark face conforming
          faceTag = zeros(2,1);
        else
          % mark faces as master and slave
          faceTag(m) = 1;
          faceTag(3-m) = 2;
        end
      end

      faceId = [faceId1,faceId2];

    end

    function faceId = getAdjacentFace(obj,index1,index2)
      % levDiff level difference between the two global index

      % faceId: local face ID of id1 cell
      % subtract the index and check in which direction
      checks = zeros(1,3);
      checks(1) = checkCellPosition(obj,index1(:,1),index2(:,1));
      checks(2) = checkCellPosition(obj,index1(:,2),index2(:,2));
      checks(3) = checkCellPosition(obj,index1(:,3),index2(:,3));

      if sum(abs(checks)) ~= 1
        % cells are not adjacent
        faceId = [];
        return
      else
        faceId = obj.mapDiffToFace([checks==1 checks==-1]);

      end



    end

    function check = checkCellPosition(obj,id1,id2)
      % check cell location along a certain direction
      % return 0: overlapping
      % return -1,1: touching
      % return 2: disjoint
      if (id1(1) <= id2(2)) && (id2(1) <= id1(2))
        check = 0;
      elseif id1(2) + 1 == id2(1)
        check = -1;
      elseif id2(2) + 1 == id1(1)
        check = 1;
      else
        check = 2;
      end

    end
  end
end

