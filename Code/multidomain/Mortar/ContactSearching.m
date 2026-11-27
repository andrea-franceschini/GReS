classdef ContactSearching < handle
  %CONTACTSEARCHING Summary of this class goes here
  %   Detailed explanation goes here

  properties
    BBTreeMaster
    BBTreeSlave
    I
    J
    cnt
    cap
    isConforming = false
    nEl
    elemConnectivity
  end

  methods
    function obj = ContactSearching(varargin)
      % CONTACTSEARCHING
      % Call as
      % ContactSearching(2DmeshSlave,2DmeshMaster,options)
      % ContactSearching(coordsSlave,coordsMaster,connSlave,connMaster,options)


      if isa(varargin{1},"Grid")
        mshSlave = varargin{MortarSide.slave};
        mshMaster = varargin{MortarSide.master};
        connMaster = mshMaster.cells.connectivity;
        connSlave =  mshSlave.cells.connectivity;
        centersMaster = mshMaster.cells.center;
        centersSlave = mshSlave.cells.center;
        coordsMaster = mshMaster.coordinates;
        coordsSlave = mshSlave.coordinates;
        k = 2;
      else
        [coordsSlave,coordsMaster,connSlave,connMaster] = deal(varargin{1:4});
        k = 4;
        centersMaster = [];
        centersSlave = [];
      end

      nSlave = size(connSlave,1);
      nMaster = size(connMaster,1);
      if nSlave == 0 || nMaster == 0
        obj.elemConnectivity = sparse(nSlave,nMaster);
        return
      end

      obj.isConforming = obj.checkConforming(coordsSlave,coordsMaster,connSlave,connMaster);
      if obj.isConforming
        return
      end

      if size(centersSlave,1) ~= nSlave || size(centersSlave,2) ~= size(coordsSlave,2)
        centersSlave = obj.computeCenters(coordsSlave,connSlave);
      end
      if size(centersMaster,1) ~= nMaster || size(centersMaster,2) ~= size(coordsMaster,2)
        centersMaster = obj.computeCenters(coordsMaster,connMaster);
      end

      obj.BBTreeMaster = BBTree(coordsMaster, connMaster, centersMaster, varargin{k+1:end});
      obj.BBTreeSlave  = BBTree(coordsSlave, connSlave, centersSlave, varargin{k+1:end});
    end


    function conn = getConnectivity(obj)
      % Rows index slave cells and columns index master cells.
      obj.findConnectivity();
      conn = obj.elemConnectivity;
    end


    function conn = getElementConnectivity(obj)
      conn = obj.getConnectivity();
    end


    function findConnectivity(obj)

      if isempty(obj.elemConnectivity) && ~isempty(obj.BBTreeSlave)

        % initial indices for sparse matrix assembly (use cap
        % for preallocation since size is unknown a priory)
        obj.cap = 1024;
        obj.I = zeros(obj.cap,1,'uint32');
        obj.J = zeros(obj.cap,1,'uint32');
        obj.cnt = 0;

        obj.tandemTraversal(1,1);

        obj.I = obj.I(1:obj.cnt);
        obj.J = obj.J(1:obj.cnt);

        % assemble the connectivity matrix
        obj.elemConnectivity = sparse(double(obj.I), double(obj.J), true, ...
          obj.BBTreeSlave.nElem, obj.BBTreeMaster.nElem);

      end

    end


    function tandemTraversal(obj,nodeSlave,nodeMaster)

      if ~obj.BBTreeSlave.intersects(nodeSlave, obj.BBTreeMaster, nodeMaster)
        return
      end

      isLeafSlave  = obj.BBTreeSlave.isLeaf(nodeSlave);
      isLeafMaster = obj.BBTreeMaster.isLeaf(nodeMaster);

      if isLeafSlave && isLeafMaster
        obj.cnt = obj.cnt + 1;

        if obj.cnt > obj.cap
          obj.cap = obj.cap * 2;
          obj.I(obj.cap,1) = 0;
          obj.J(obj.cap,1) = 0;
        end

        obj.I(obj.cnt) = obj.BBTreeSlave.leafElem(nodeSlave);
        obj.J(obj.cnt) = obj.BBTreeMaster.leafElem(nodeMaster);
        return
      end

      if ~isLeafSlave && ~isLeafMaster
        cS = obj.BBTreeSlave.children(nodeSlave,:);
        cM = obj.BBTreeMaster.children(nodeMaster,:);

        tandemTraversal(obj,cS(1), cM(1));
        tandemTraversal(obj,cS(1), cM(2));
        tandemTraversal(obj,cS(2), cM(1));
        tandemTraversal(obj,cS(2), cM(2));

      elseif ~isLeafSlave
        cS = obj.BBTreeSlave.children(nodeSlave,:);
        tandemTraversal(obj,cS(1), nodeMaster);
        tandemTraversal(obj,cS(2), nodeMaster);

      else
        cM = obj.BBTreeMaster.children(nodeMaster,:);
        tandemTraversal(obj,nodeSlave, cM(1));
        tandemTraversal(obj,nodeSlave, cM(2));
      end
    end


    function isConforming = checkConforming(obj,X1,X2,T1,T2)
      isConforming = false;
      if ~isequal(size(X1),size(X2)) || size(T1,1) ~= size(T2,1)
        return
      end

      [X1s,p1] = sortrows(X1);
      [X2s,p2] = sortrows(X2);
      if ~isequal(X1s,X2s)
        return
      end

      [T1s,e1] = obj.canonicalCells(T1,p1);
      [T2s,e2] = obj.canonicalCells(T2,p2);
      isConforming = isequal(T1s,T2s);
      if isConforming
        n = size(T1,1);
        obj.elemConnectivity = sparse(e1,e2,true,n,n);
      end
    end

  end


  methods (Static, Access = private)

    function [cells,order] = canonicalCells(conn,nodeOrder)
      if isa(conn,'ArrayOfArrays')
        conn = conn.toMatrix();
      end
      nodeIds = zeros(numel(nodeOrder),1);
      nodeIds(nodeOrder) = 1:numel(nodeOrder);
      ids = zeros(size(conn));
      active = conn > 0;
      ids(active) = nodeIds(conn(active));
      ids = ids(:,any(active,1));
      [cells,order] = sortrows(sort(ids,2));
    end

  end


  methods (Static)


    function centers = computeCenters(coordinates, connectivity)

      if isa(connectivity, 'ArrayOfArrays')
        v = connectivity.getData();
        nVert = connectivity.arraySize();
      else
        poly = connectivity';
        v = poly(poly > 0);
        nVert = sum(connectivity > 0, 2);
      end

      poly = coordinates(v,:);
      centers = computePolygonCentroid(poly,nVert);

    end


  end


end

