classdef BBTree < handle
  %BBTREE Bounding-volume binary tree for polytopal contact search.
  %
  % General-purpose tree built from:
  %   coordinates   : [nNodes x dimEmbed]
  %   connectivity  : [nElem x nNodePerElem] numeric matrix
  %                   or ArrayOfArrays
  %   centers       : [nElem x dim]
  %
  % Tree nodes store k-DOP / k-top primitive intervals:
  %   nodeBoxes(i,:) = [min_1 ... min_k max_1 ... max_k]
  %
  % Leaves correspond to single elements and are identified by leafElem(i).

  properties (SetAccess = private)
    coordinates
    connectivity
    centers
    nElem
    dim
    polytop
    scale = 0.025

    children        % [nNodesTree x 2], zero if leaf
    nodeBoxes       % [nNodesTree x 2*nPrim]
    leafElem        % [nNodesTree x 1], zero for internal nodes
  end

  methods
    function obj = BBTree(coordinates, connectivity, centers, varargin)
      % BBTree(coordinates, connectivity, centers)
      % BBTree(..., 'Scale', val, 'Polytop', P)

      obj.coordinates  = coordinates;
      obj.connectivity = connectivity;
      obj.centers      = centers;

      obj.nElem = size(centers,1);
      obj.dim   = size(centers,2);

      % defaults
      default = struct('scale',0.025,'polytop',BBTree.defaultPolytop(obj.dim));
      parm = readInput(default,varargin{:});

      obj.scale = parm.scale;
      obj.polytop = parm.polytop;

      obj.build();
    end




    function tf = isLeaf(obj, nodeId)
      tf = obj.children(nodeId,1) == 0;
    end



    function tf = intersects(obj, nodeId1, otherTree, nodeId2)
      % Check node-node intersection between this tree and another tree

      box1 = obj.nodeBoxes(nodeId1,:);
      box2 = otherTree.nodeBoxes(nodeId2,:);

      nPrim = size(obj.polytop,2);

      min1 = box1(1:nPrim);
      max1 = box1(nPrim+1:end);
      min2 = box2(1:nPrim);
      max2 = box2(nPrim+1:end);

      tf = ~any(min1 > max2 | min2 > max1);
    end
  end



  methods (Access = private)
    function build(obj)
      % Build tree iteratively using a queue of element subsets.
      %

      nTreeMax = max(1, 2*obj.nElem - 1);
      nPrim    = size(obj.polytop,2);

      obj.children = zeros(nTreeMax, 2);
      obj.nodeBoxes = zeros(nTreeMax, 2*nPrim);
      obj.leafElem = zeros(nTreeMax, 1);

      % queue of element lists per node
      elemSets = cell(nTreeMax, 1);
      elemSets{1} = (1:obj.nElem).';

      nextFree = 2;
      nodeId   = 1;

      while nodeId < nextFree
        elemId = elemSets{nodeId};

        [box, leftElem, rightElem] = obj.populateNode(elemId);
        obj.nodeBoxes(nodeId,:) = box;

        if isempty(leftElem)
          % leaf
          obj.leafElem(nodeId) = elemId;
        else
          lId = nextFree;
          rId = nextFree + 1;

          obj.children(nodeId,:) = [lId, rId];
          elemSets{lId} = leftElem;
          elemSets{rId} = rightElem;

          nextFree = nextFree + 2;
        end

        nodeId = nodeId + 1;
      end

      % trim unused preallocation
      lastUsed = nextFree - 1;
      obj.children = obj.children(1:lastUsed,:);
      obj.nodeBoxes = obj.nodeBoxes(1:lastUsed,:);
      obj.leafElem = obj.leafElem(1:lastUsed);
    end



    function [box, leftElem, rightElem] = populateNode(obj, elemId)
      % Build node bounding box and split element subset if needed.

      connSel = getRows(obj.connectivity, elemId);

      if isa(connSel, 'ArrayOfArrays')
        [flatConn, ~] = connSel.getData();
      else
        flatConn = connSel.';
        flatConn = flatConn(:);
        flatConn = flatConn(flatConn > 0);
      end

      nodes = unique(flatConn);
      coords = obj.coordinates(nodes, 1:obj.dim);

      prim = coords * obj.polytop;
      pmin = min(prim, [], 1);
      pmax = max(prim, [], 1);

      delta = max(pmax - pmin);
      if delta > 0
        pmin = pmin - obj.scale * delta;
        pmax = pmax + obj.scale * delta;
      end

      box = [pmin, pmax];

      nElemLoc = numel(elemId);

      if nElemLoc == 1
        leftElem  = [];
        rightElem = [];
        return
      end

      if nElemLoc == 2
        leftElem  = elemId(1);
        rightElem = elemId(2);
        return
      end

      % Split along widest primitive direction
      [~, splitDir] = max(pmax - pmin);

      cproj = obj.centers(elemId,1:obj.dim) * obj.polytop(:,splitDir);

      pivot = median(cproj);
      goLeft = cproj < pivot;

      leftElem  = elemId(goLeft);
      rightElem = elemId(~goLeft);

      % fallback: mean
      if isempty(leftElem) || isempty(rightElem)
        pivot = mean(cproj);
        goLeft = cproj < pivot;
        leftElem  = elemId(goLeft);
        rightElem = elemId(~goLeft);
      end

      % last fallback: deterministic half split
      if isempty(leftElem) || isempty(rightElem)
        n = nElemLoc;
        ord = 1:n;
        mid = ceil(n/2);
        leftElem  = elemId(ord(1:mid));
        rightElem = elemId(ord(mid+1:end));
      end
    end
  end




  methods (Static, Access = private)
    function P = defaultPolytop(dim)
      switch dim
        case 2
          % 8-DOP written as dim x nPrim
          P = [ 1  0 -1  1;
            0  1  1  1 ];
        case 3
          % 18-DOP written as dim x nPrim
          P = [ 1  0  0  1  1  0  1  1  0;
            0  1  0  1  0  1 -1  0  1;
            0  0  1  0  1  1  0 -1 -1 ];
        otherwise
          error('BBTree:UnsupportedDimension', ...
            'Only dim = 2 or dim = 3 are supported.');
      end
    end
  end
end

