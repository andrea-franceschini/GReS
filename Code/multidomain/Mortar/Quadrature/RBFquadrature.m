classdef RBFquadrature < MortarQuadrature
  % Radial Basis Function class
  % Used to interpolate master basis functions into the slave side

  properties
    activeGPmap        % index map to access gp related information
    nInt               % number of interpolation points per master element
    nGP                % number of gauss points in input
    rbfType
    detJw
  end

  properties (Access = private)
    wF
    w1
    ptsRBF

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
    function obj = RBFquadrature(multType,grids,input)
      %
      obj@MortarQuadrature(multType,grids,input);
      default = struct('nInt',5,'RBFtype',"gauss");
      params = readInput(default,input);
      obj.nInt = params.nInt;
      obj.rbfType = params.RBFtype;
      obj.getWeights();
    end

  end

  methods (Access = public)

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


    function isPairActive = processMortarPair(obj,is,im,k)

      isPairActive = true;

      posGP = obj.gpsCoord;

      elem = getElem(obj,1,im);

      if class(elem)=="Triangle"
        nIntPts = sum(1:obj.nInt);
      else
        nIntPts = obj.nInt^2;
      end

      ptsInt = obj.ptsRBF(1:nIntPts,[3*im-2 3*im-1 3*im]);

      [fiNM,id1] = obj.computeRBFfiNM(ptsInt,posGP);
      % id1: flags gp that are too far from master support

      % compute interpolated gp coordinates
      xiMaster = (fiNM*obj.wF(:,[2*im-1 2*im]))./(fiNM*obj.w1(:,im));

      obj.suppFlag = elem.checkInRange(xiMaster,1e-3);

      if ~any(obj.suppFlag)
        isPairActive = false;
        return
      end

      % store infos in maps
      obj.interfacePairs(k,:) = [is im];

      nGsupp = sum(obj.suppFlag);
      gpId = (is-1)*obj.maxGP + obj.countGP;
      obj.gpCoords{1}(gpId+1:gpId+nGsupp,1) = xiMaster(obj.suppFlag,1);
      obj.gpCoords{1}(gpId+1:gpId+nGsupp,2) = xiMaster(obj.suppFlag,2);
      obj.gpCoords{2}(gpId+1:gpId+nGsupp,1) = obj.gpsCoordLoc(obj.suppFlag,1);
      obj.gpCoords{2}(gpId+1:gpId+nGsupp,2) = obj.gpsCoordLoc(obj.suppFlag,2);
      obj.detJw(gpId+1:gpId+nGsupp) = obj.dJwSlave(obj.suppFlag);

    end

    function finalizeMortarMaps(obj)
      % remove useless entries after mortar preallocation
      id = ~any(obj.interfacePairs,2);
      id2 = obj.detJw == 0;

      obj.activeGPmap([false;id]) = [];
      obj.interfacePairs(id,:) = [];

      obj.detJw(id2) = [];
      for i = 1:2
        obj.gpCoords{i}(id2,:) = [];
      end

      obj.numbInterfacePairs = size(obj.interfacePairs,1);
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

  end

  methods (Access = private)
    %
    function getWeights(obj)

      m = MortarSide.master;
      surf = obj.grids(m).surfaces;
      coordinates = obj.grids(m).coordinates;

      vtk = surf.VTKType;
      nTri = sum(vtk == VTKType.Tri);
      nQuad = sum(any([vtk == VTKType.Quad; vtk == VTKType.Quad9]));
      numPts =  nTri*sum(1:obj.nInt) + nQuad*(obj.nInt)^2;

      nElM = surf.num;
      weighF = zeros(numPts,2*nElM);
      weigh1 = zeros(numPts,nElM);
      pts = zeros(numPts,nElM*3);

      % solve local system to get weight of interpolant
      warning('off','MATLAB:nearlySingularMatrix')

      k = 0;

      for vtkMaster = obj.grids(m).surfaces.vtkTypes

        elem = FiniteElementType.create(vtkMaster,...
          obj.grids(m),'gaussOrder',obj.gaussOrder);
        listMaster = obj.grids(m).getSurfByVTKId(vtkMaster);

        topol = obj.grids(m).getCellNodes(listMaster);


        for im = listMaster'

          nodes = topol(im,:);
          coord = coordinates(nodes,:);

          f = getInterpolationPoints(obj,elem);
          N = elem.computeBasisF(f);
          % get coords of interpolation points in the real space
          pos = N*coord;


          nptInt = size(pos,1);

          fiMM = obj.computeRBFfiMM(pos);

          x = fiMM\[f ones(size(ptsInt,1),1)];

          weighF(1:nptInt,k+1:k+2) = x(:,1:2);
          weigh1(1:nptInt,im) = x(:,3);

          pts(1:nptInt,[3*im-2 3*im-1 3*im]) = ptsInt;
          k = k+2;
        end

      end

      obj.wF = weighF;
      obj.w1 = weigh1;
      obj.ptsRBF = pts;

    end

    function fiMM = computeRBFfiMM(obj,ptsInt)
      r = obj.computeRBFradius(ptsInt);
      d = sqrt((ptsInt(:,1) - ptsInt(:,1)').^2 + (ptsInt(:,2) - ptsInt(:,2)').^2 + (ptsInt(:,3) - ptsInt(:,3)').^2);
      fiMM = obj.rbfInterp(d,r,obj.rbfType);
    end

    function [fiNM,id] = computeRBFfiNM(obj,ptsInt,ptsGauss)
      % id: id of points that has a distance < r with at least one
      % master points (to autmatically sort out gauss points that are too
      % far from the master interpolation points)
      d = sqrt((ptsGauss(:,1) - ptsInt(:,1)').^2 + (ptsGauss(:,2) - ptsInt(:,2)').^2 + (ptsGauss(:,3) - ptsInt(:,3)').^2);
      r = sqrt((max(ptsInt(:,1)) - min(ptsInt(:,1)))^2 + (max(ptsInt(:,2)) - min(ptsInt(:,2)))^2 + (max(ptsInt(:,3)) - min(ptsInt(:,3)))^2);
      id = ~all(d>=r,2);
      fiNM = obj.rbfInterp(d,r,obj.rbfType);
    end

    % function [intPts,pos] = getRBFfunction(obj,elId,elem)
    %   % evaluate shape function in the real space and return position of
    %   % integration points in the real space
    %   surf = obj.grids(MortarSide.master).surface;
    %   surfNodes = msh.surfaces(elId,:);
    %   coord = msh.coordinates(surfNodes,:);
    %   elem = obj.getElem(1,elId);
    %   % place interpolation points in a regular grid
    %   intPts = getInterpolationPoints(obj,elem);
    %   bf = elem.computeBasisF(intPts);
    %   % get coords of interpolation points in the real space
    %   pos = bf*coord;
    % end

    function intPts = getInterpolationPoints(obj,vtk)
      switch vtk
        case VTKType.Tri
          % uniform grid in triangle
          intPts = zeros(sum(1:obj.nInt),2);
          p = linspace(0,1,obj.nInt);
          k = obj.nInt;
          c = 0;
          for i = 1:obj.nInt
            intPts(c+1:c+k,1) = p(1:k)';
            c = c + k;
            k = k-1;
          end
          intPts(:,2) = repelem(p,obj.nInt:-1:1);
        case {VTKType.Quad,VTKType.Quad9}
          intPts = linspace(-1,1, obj.nInt);
          [y, x] = meshgrid(intPts, intPts);
          intPts = [x(:), y(:)];
      end
    end
  end


  methods (Static)

    function rbf = rbfInterp(d,r,type)
      % compute row of rbf interpolation matrix
      switch type
        case 'gauss'
          rbf = exp(-d.^2/r^2);
        case 'wendland'
          v = 1-d./r;
          v(v<0) = 0;
          rbf = v.^4.*(1+d./r);
      end
    end

    function r = computeRBFradius(ptsInt)
      % compute the radius of the RBF interpolation based on
      % coordinates of interpolation points
      r = sqrt((max(ptsInt(:,1)) - min(ptsInt(:,1)))^2 + ...
        (max(ptsInt(:,2)) - min(ptsInt(:,2)))^2 + ...
        (max(ptsInt(:,3)) - min(ptsInt(:,3)))^2);
    end
  end

end

