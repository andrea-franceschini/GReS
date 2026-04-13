classdef ContactSearching
  %CONTACTSEARCHobj.ING Summary of this class goes here
  %   Detailed explanation goes here

  properties
    BBTreeMaster
    BBTreeSlave
    I
    J
    cnt
    cap
  end

  methods
    function obj = ContactSearching(varargin)
      % CONTACTSEARCHING
      % Call as
      % ContactSearching(2DmeshMaster,2DmeshSlave,options)
      % ContactSearching(topolMaster,topolSlave,options)


      if isMatrix(varargin{1})
        mshMaster = varargin{2};
        mshSlave = varargin{3};
        connMaster = mshMaster.surfaces.connectivity;
        connSlave =  mshSlave.surfaces.connectivity;
        centersMaster = mshMaster.surfaces.center;
        centersSlave = mshSlave.surfaces.center;
        coordsMaster = mshMaster.coordinates;
        coordsSlave = mshSlave.coordinates;
        k = 2;
      else
        [coordsMaster,coordsSlave,connMaster,connSlave] = deal(varargin{1:4});
        k = 4;
        centersMaster = computeSurfaceCenters(coordsMaster, connMaster);
        centersSlave  = computeSurfaceCenters(coordsSlave,  connSlave);
      end

      opts = varargin{k+1:end};

      % build bounding-box trees
      obj.BBTreeMaster = BBTree(coordsMaster, connMaster, centersMaster, opts);
      obj.BBTreeSlave  = BBTree(coordsSlave, connSlave, centersSlave, opts);

    end


    function elemConnectivity = getElementConnectivity(obj)


      % initial indices for sparse matrix assembly (use cap
      % for preallocation since size is unknown a priory)
      obj.cap = 1024;
      obj.I = zeros(obj.cap,1,'uint32');
      obj.J = zeros(obj.cap,1,'uint32');
      obj.cnt = 0;

      tandemTraversal(1,1);

      obj.I = obj.I(1:obj.cnt);
      obj.J = obj.J(1:obj.cnt);

      % assemble the connectivity matrix
      elemConnectivity = sparse(double(obj.I), double(obj.J), true, ...
        obj.BBTreeSlave.nElem, obj.BBTreeMaster.nElem);


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

        tandemTraversal(cS(1), cM(1));
        tandemTraversal(cS(1), cM(2));
        tandemTraversal(cS(2), cM(1));
        tandemTraversal(cS(2), cM(2));

      elseif ~isLeafSlave
        cS = obj.BBTreeSlave.children(nodeSlave,:);
        tandemTraversal(cS(1), nodeMaster);
        tandemTraversal(cS(2), nodeMaster);

      else
        cM = obj.BBTreeMaster.children(nodeMaster,:);
        tandemTraversal(nodeSlave, cM(1));
        tandemTraversal(nodeSlave, cM(2));
      end
    end


    function centers = computeSurfaceCenters(coordinates, connectivity)

      if isa(connectivity, 'ArrayOfArrays')
        v = connectivity.getData();
        nVert = connectivity.arraySize();
      else
        poly = connectivity';
        v = poly(:);
        nVert = sum(connectivity > 0, 2);
      end

      poly = coordinates(v,:);
      centers = computePolygonCentroid(poly,nVert);

    end

  end
end

