classdef ElementBasedQuadrature < MortarQuadrature
  % Implement utilities to perform element based integration given a pair
  % of master and slave element in contact
  % Can be equipped with the RBF scheme to cheaply sort out elements that
  % are not really connected

  % REFS: Puso,2004, A mortar segment-to-segment contact method for large
  % deformation solid mechanics

  properties
    activeGPmap        % index map to access gp related information
    nGP                % number of gauss points in input
    detJw
  end

  properties (Access = private)
    % element based infos
    idSlave          % current slave element being processed
    gpsCoord     % current list of 3D gp coordinates for RBF interpolation
    gpsCoordLoc  % current list of local gp coordinates
    dJwSlave
    suppFlag
    countGP
    maxGP
  end

  methods
    function obj = ElementBasedQuadrature(mortar,nGP)
      obj@MortarQuadrature(mortar,nGP);
      obj.nGP = nGP;
    end


    function processMortarPairs(obj)

      % initialize the maps to store mortar quadrature info
      obj.maxGP = obj.nGP^2;
      %       nConnections = nnz(obj.mortar.mesh.elemConnectivity);
      totGP = obj.msh(2).nSurfaces*obj.maxGP;
      obj.gpCoords = {zeros(totGP,2);
        zeros(totGP,2)};

      nConnections = nnz(obj.mortar.mesh.elemConnectivity);

      obj.mortarPairs = zeros(nConnections,2);
      obj.detJw = zeros(totGP,1);

      obj.activeGPmap = zeros(nConnections+1,1);

      nM = full(sum(obj.mortar.mesh.elemConnectivity,1));
      nM = [0 cumsum(nM)];

      for is = 1:obj.msh(2).nSurfaces

        % reset slave element based info
        elemSlave = obj.mortar.getElem(2,is);
        obj.idSlave = is;
        obj.countGP = 0;
        obj.gpsCoord = getGPointsLocation(elemSlave,is);
        obj.gpsCoordLoc = elemSlave.GaussPts.coord;
        obj.dJwSlave = getDerBasisFAndDet(elemSlave,is);

        imList = find(obj.mortar.mesh.elemConnectivity(:,is));

        for j = 1:numel(imList)
          im = imList(j);
          k = nM(is)+ j;
          % k (global index to write without race conditions)
          isPairActive = processMortarPair(obj,is,im,k);
          if ~isPairActive
            obj.activeGPmap(k+1) = obj.activeGPmap(k);
            continue
          end

          % sort out gauss points already projected
          obj.activeGPmap(k+1) = obj.activeGPmap(k) + sum(obj.suppFlag);
          obj.gpsCoord = obj.gpsCoord(~obj.suppFlag,:);
          obj.gpsCoordLoc = obj.gpsCoordLoc(~obj.suppFlag,:);
          obj.dJwSlave = obj.dJwSlave(~obj.suppFlag);
          obj.countGP = obj.countGP + sum(obj.suppFlag);


        end

        if obj.countGP ~= elemSlave.GaussPts.nNode
          warning("Some gauss point not projected for element %i",is)
        end

      end

      finalizeMortarMaps(obj);

    end



    function isPairActive = processMortarPair(obj,is,im,k)

      isPairActive = true;

      % compute interpolated gp coordinates and update supportFlag
      xiMaster = projectGP(obj,is,im);


      if ~any(obj.suppFlag)
        isPairActive = false;
        return
      end

      % store infos in maps
      obj.mortarPairs(k,:) = [is im];

      nGsupp = sum(obj.suppFlag);
      gpId = (is-1)*obj.maxGP + obj.countGP;
      obj.gpCoords{1}(gpId+1:gpId+nGsupp,1) = xiMaster(obj.suppFlag,1);
      obj.gpCoords{1}(gpId+1:gpId+nGsupp,2) = xiMaster(obj.suppFlag,2);
      obj.gpCoords{2}(gpId+1:gpId+nGsupp,1) = obj.gpsCoordLoc(obj.suppFlag,1);
      obj.gpCoords{2}(gpId+1:gpId+nGsupp,2) = obj.gpsCoordLoc(obj.suppFlag,2);
      obj.detJw(gpId+1:gpId+nGsupp) = obj.dJwSlave(obj.suppFlag);

    end

    function xiM = projectGP(obj,is,im)
      % xi: reference coordinates of the gauss point
      % get nodal normal
      elM = getElem(obj.mortar,1,im);
      elS = getElem(obj.mortar,2,is);
      nodeS = obj.msh(2).surfaces(is,:);
      X = obj.gpsCoord;                    % real position of gauss pts
      xiS = obj.gpsCoordLoc;

      ngp = size(xiS,1);

      obj.suppFlag = false(ngp,1);

      xiM = zeros(ngp,2);
      itMax = 8;
      tol = 1e-10;

      for i = 1:ngp
        Ns = elS.computeBasisF(xiS(i,:));
        ng = Ns*obj.mortar.mesh.avgNodNormal{2}(nodeS,:); % slave normal at GP
        ng = ng/norm(ng);
        xiM(i,:) = elS.centroid;                          % initial guess
        iter = 0;
        w = 0;
        nodeM = obj.msh(1).surfaces(im,:);
        coordM = obj.msh(1).coordinates(nodeM,:);
        Nm = elM.computeBasisF(xiM(i,:));
        rhs = (Nm*coordM - w*ng - X(i,:))';
        %
        while (norm(rhs,2) > tol) && (iter < itMax)
          iter = iter+1;
          dN = elM.computeDerBasisF(xiM(i,:));
          J1 = dN*coordM;
          J = [J1' -ng'];
          ds = J\(-rhs);
          xiM(i,:) = xiM(i,:) + (ds(1:2))';
          w = w + ds(3);
          Nm =  elM.computeBasisF(xiM(i,:));
          rhs = (Nm*coordM - w*ng - X(i,:))';
        end

        % check if gp lies in master reference space and update gp flag
        if elM.checkInRange(xiM(i,:))
          obj.suppFlag(i) = true;
        end
      end
    end


    function gpCoords = getSlaveGPCoords(obj,kPair)
      i1 = obj.activeGPmap(kPair);
      i2 = obj.activeGPmap(kPair+1);
      gpCoords = obj.gpCoords{2}(i1+1:i2,:);
    end

    function gpCoords = getMasterGPCoords(obj,kPair)
      i1 = obj.activeGPmap(kPair);
      i2 = obj.activeGPmap(kPair+1);
      gpCoords = obj.gpCoords{1}(i1+1:i2,:);
    end


    function dJweighed = getIntegrationWeights(obj,kPair)
      i1 = obj.activeGPmap(kPair);
      i2 = obj.activeGPmap(kPair+1);
      dJweighed = obj.detJw(i1+1:i2);
    end


    function finalizeMortarMaps(obj)
      % remove useless entries after mortar preallocation
      id = ~any(obj.mortarPairs,2);
      id2 = obj.detJw == 0;

      obj.activeGPmap([false;id]) = [];
      obj.mortarPairs(id,:) = [];

      obj.detJw(id2) = [];
      for i = 1:2
        obj.gpCoords{i}(id2,:) = [];
      end

      obj.numbMortarPairs = size(obj.mortarPairs,1);
    end


  end
end

