classdef ElementBasedQuadrature < MortarQuadrature
  % Implement utilities to perform element based integration given a pair
  % of master and slave element in contact
  % Can be equipped with the RBF scheme to cheaply sort out elements that
  % are not really connected

  % REFS: Puso,2004, A mortar segment-to-segment contact method for large
  % deformation solid mechanics

  properties
    activeGPmap        % index map to access gp related information
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
    function obj = ElementBasedQuadrature(multType,grids,input)

      obj@MortarQuadrature(multType,grids,input);

    end


    function processMortarPairs(obj,connectivity)

      % initialize the maps to store mortar quadrature info
      obj.maxGP = Gauss.getNumPtsFromOrder(obj.gaussOrder)^2;

      s = MortarSide.slave;
      m = MortarSide.master;
      totGP = obj.grids(s).surfaces.num*obj.maxGP;
      obj.gpCoords = {zeros(totGP,2);
                      zeros(totGP,2)};

      nConnections = nnz(connectivity);

      obj.interfacePairs = zeros(nConnections,2);
      obj.detJw = zeros(totGP,1);

      obj.activeGPmap = zeros(nConnections+1,1);


      nM = full(sum(connectivity,2));
      nM = [0;cumsum(nM)];




      for vtkSlave = obj.grids(s).surfaces.vtkTypes

        elemSlave = FiniteElementType.create(vtkSlave,obj.grids(s));
        listSlave = obj.grids(s).getSurfByVTKId(vtkSlave);

        for vtkMaster = obj.grids(m).surfaces.vtkTypes

          elemMaster = FiniteElementType.create(vtkMaster,obj.grids(m));
          listMaster = obj.grids(m).getSurfByVTKId(vtkMaster);

          mask = false(obj.grids(m).surfaces.num,1);
          mask(listMaster) = true;

          for is = listSlave'

            obj.idSlave = is;
            obj.countGP = 0;
            obj.gpsCoord = getGPointsLocation(elemSlave,is);
            obj.gpsCoordLoc = elemSlave.getGauss.coord;
            obj.dJwSlave = getDerBasisFAndDet(elemSlave,is);


            % get master elements belonging to that vtk type
            imList = find(connectivity(is,:));
            imList = imList(mask(imList));

            for j = 1:numel(imList)
              im = imList(j);
              k = nM(is)+ j;
              % k (global index to write without race conditions)
              isPairActive = processMortarPair(obj,is,im,elemSlave,elemMaster,k);

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

            if obj.countGP ~= elemSlave.getGauss.nNode
              warning("Some gauss point not projected for element %i",is)
            end

          end
        end
      end


      finalizeMortarMaps(obj);
      computeAreaSlave(obj);

    end



    function isPairActive = processMortarPair(obj,is,im,elemSlave,elemMaster,k)

      isPairActive = true;

      s = MortarSide.slave;
      m = MortarSide.master;

      % compute interpolated gp coordinates and update supportFlag
      xiMaster = projectGP(obj,is,im,elemSlave,elemMaster);


      if ~any(obj.suppFlag)
        isPairActive = false;
        return
      end

      % store infos in maps
      obj.interfacePairs(k,[s m]) = [is im];

      nGsupp = sum(obj.suppFlag);
      gpId = (is-1)*obj.maxGP + obj.countGP;
      obj.gpCoords{m}(gpId+1:gpId+nGsupp,1) = xiMaster(obj.suppFlag,1);
      obj.gpCoords{m}(gpId+1:gpId+nGsupp,2) = xiMaster(obj.suppFlag,2);
      obj.gpCoords{s}(gpId+1:gpId+nGsupp,1) = obj.gpsCoordLoc(obj.suppFlag,1);
      obj.gpCoords{s}(gpId+1:gpId+nGsupp,2) = obj.gpsCoordLoc(obj.suppFlag,2);
      obj.detJw(gpId+1:gpId+nGsupp) = obj.dJwSlave(obj.suppFlag);

    end

    function xiM = projectGP(obj,is,im,elS,elM)
      % xi: reference coordinates of the gauss point
      % get nodal normal

      gridSlave = obj.grids(MortarSide.slave);
      gridMaster = obj.grids(MortarSide.master);
      nodeS = gridSlave.getSurfNodes(is);
      nodeM =  gridMaster.getSurfNodes(im);
      coordM = gridMaster.coordinates(nodeM,:);
      X = obj.gpsCoord;                    % real position of gauss pts
      xiS = obj.gpsCoordLoc;

      normal = gridSlave.surfaces.avgNodNormal(nodeS,:);

      ngp = size(xiS,1);

      obj.suppFlag = false(ngp,1);

      xiM = zeros(ngp,2);
      itMax = 8;
      tol = 1e-10;

      for i = 1:ngp
        Ns = elS.computeBasisF(xiS(i,:));
        ng = Ns*normal;        % slave normal at GP
        ng = ng/norm(ng);
        xiM(i,:) = elS.centroid;                        % initial guess for gp
        iter = 0;
        w = 0;
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
      gpCoords = obj.gpCoords{MortarSide.slave}(i1+1:i2,:);
    end

    function gpCoords = getMasterGPCoords(obj,kPair)
      i1 = obj.activeGPmap(kPair);
      i2 = obj.activeGPmap(kPair+1);
      gpCoords = obj.gpCoords{MortarSide.master}(i1+1:i2,:);
    end


    function dJweighed = getIntegrationWeights(obj,kPair)
      i1 = obj.activeGPmap(kPair);
      i2 = obj.activeGPmap(kPair+1);
      dJweighed = obj.detJw(i1+1:i2);
    end


    function finalizeMortarMaps(obj)
      % remove useless entries after mortar preallocation
      id = ~any(obj.interfacePairs,2);
      id2 = obj.detJw == 0;

      obj.activeGPmap([false;id]) = [];
      obj.interfacePairs(id,:) = [];

      obj.detJw(id2) = [];
      for side = [MortarSide.master,MortarSide.slave]
        obj.gpCoords{side}(id2,:) = [];
      end

      obj.numbInterfacePairs = size(obj.interfacePairs,1);
    end


  end
end

