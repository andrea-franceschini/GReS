function elemConnectivity = contactSearch( ...
  coordMaster, connMaster, ...
  coordSlave,  connSlave, varargin)
%CONTACTSEARCH Contact search between two surface meshes using BB-trees.
%
% elemConnectivity = contactSearch( ...
%     coordinatesMaster, connectivityMaster, ...
%     coordinatesSlave,  connectivitySlave)
%
% Inputs
%   coordinatesMaster : [nNodesM x dim]
%   connectivityMaster: [nElemM x nNodePerElem] matrix or ArrayOfArrays
%   coordinatesSlave  : [nNodesS x dim]
%   connectivitySlave : [nElemS x nNodePerElem] matrix or ArrayOfArrays
%
% Optional name-value pairs
%   'Scale'   : bounding volume expansion factor (default 0.025)
%   'Polytop' : custom primitive directions matrix [dim x nPrim]
%
% Output
%   elemConnectivity  : sparse(nElemSlave, nElemMaster)
%
% Convention
%   Rows    -> slave elements
%   Columns -> master elements

centersMaster = computeSurfaceCenters(coordMaster, connMaster);
centersSlave  = computeSurfaceCenters(coordSlave,  connSlave);

% build bounding-box trees
treeMaster = BBTree(coordMaster, connMaster, centersMaster, varargin{:});
treeSlave  = BBTree(coordSlave, connSlave, centersSlave, varargin{:});

% initial indices for sparse matrix assembly (use cap 
% for preallocation since size is unknown a priory)
cap = 1024;
I = zeros(cap,1,'uint32');
J = zeros(cap,1,'uint32');
cnt = 0;

tandemTraversal(1,1);

I = I(1:cnt);
J = J(1:cnt);

% assemble the connectivity matrix
elemConnectivity = sparse(double(I), double(J), true, ...
    treeSlave.nElem, treeMaster.nElem);



% nested function for tandem traversal
function tandemTraversal(nodeSlave, nodeMaster)

  if ~treeSlave.intersects(nodeSlave, treeMaster, nodeMaster)
    return
  end

  isLeafSlave  = treeSlave.isLeaf(nodeSlave);
  isLeafMaster = treeMaster.isLeaf(nodeMaster);

  if isLeafSlave && isLeafMaster
    cnt = cnt + 1;

    if cnt > cap
      cap = cap * 2;
      I(cap,1) = 0;
      J(cap,1) = 0;
    end

    I(cnt) = treeSlave.leafElem(nodeSlave);
    J(cnt) = treeMaster.leafElem(nodeMaster);
    return
  end

  if ~isLeafSlave && ~isLeafMaster
    cS = treeSlave.children(nodeSlave,:);
    cM = treeMaster.children(nodeMaster,:);

    tandemTraversal(cS(1), cM(1));
    tandemTraversal(cS(1), cM(2));
    tandemTraversal(cS(2), cM(1));
    tandemTraversal(cS(2), cM(2));

  elseif ~isLeafSlave
    cS = treeSlave.children(nodeSlave,:);
    tandemTraversal(cS(1), nodeMaster);
    tandemTraversal(cS(2), nodeMaster);

  else
    cM = treeMaster.children(nodeMaster,:);
    tandemTraversal(nodeSlave, cM(1));
    tandemTraversal(nodeSlave, cM(2));
  end
end

end




% helper function for computing the surface centroid
