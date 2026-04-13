classdef SegmentBasedQuadrature < MortarQuadrature
  % Implement utilities to perform segment based integration given a pair
  % of master and slave element in contact
  % Can be equipped with the RBF scheme to cheaply sort out elements that
  % are not really connected

  % REFS: Puso,2004, A mortar segment-to-segment contact method for large
  % deformation solid mechanics
  
  properties
    detJtri
    triGaussOrder
  end

  properties (Access = private)
    wTri      % store weight for gauss integration on triangles
    maxTriPerPair
  end
  
  methods
    function obj = SegmentBasedQuadrature(multType,grids,input)

      obj@MortarQuadrature(multType,grids);

      input = readInput(struct('gaussOrder',3),input);
      obj.triGaussOrder = input.gaussOrder;

      gaussTri = Gauss(5,obj.triGaussOrder);      % 5 is the vtk type of triangles
      obj.wTri = gaussTri.weight;

      isQuadratic = [any(grids(1).surfaces.VTKType == 28);
                     any(grids(2).surfaces.VTKType == 28)];

      if sum(isQuadratic)==0 % no quadratic elements
        obj.maxTriPerPair = 6;
      elseif sum(isQuadratic)==1 % quadratic on one side
        obj.maxTriPerPair = 24;
      elseif sum(isQuadratic)==2 % quadratic on both sides
        obj.maxTriPerPair = 96;
      end


    end

    function processMortarPairs(obj,connectivity)

      % initialize the maps to store mortar quadrature infos
      nConnections = nnz(connectivity);
      totTri = nConnections*obj.maxTriPerPair;
      ng = numel(obj.wTri);
      obj.gpCoords = {zeros(totTri,ng,2);
        zeros(totTri,ng,2)};
      obj.interfacePairs = zeros(totTri,2);
      obj.detJtri = zeros(totTri,1);

      nM = full(sum(connectivity,2));
      nM = [0;cumsum(nM)];

      gridS = obj.grids(MortarSide.slave);
      gridM = obj.grids(MortarSide.master);

      % process homogenous vtk types on both sides of the grid

      for vtkSlave = gridS.surfaces.vtkTypes

        elemSlave = FiniteElementType.create(vtkSlave,gridS);
        listSlave = gridS.getSurfByVTKId(vtkSlave);

        for vtkMaster = gridM.surfaces.vtkTypes

          elemMaster = FiniteElementType.create(vtkMaster,gridM);
          listMaster = gridM.getSurfByVTKId(vtkMaster);

          mask = false(gridM.surfaces.num,1);
          mask(listMaster) = true;

          for is = listSlave'

            % get master elements belonging to that vtk type
            imList = find(connectivity(is,:));
            imList = imList(mask(imList));

            for j = 1:numel(imList)
              im = imList(j);
              k = nM(is)+ j;
              % k (global index to write without race conditions)
              processMortarPair(obj,is,im,elemSlave,elemMaster,k);
            end
          end
        end
      end

      finalizeMortarMaps(obj);
      computeAreaSlave(obj);

    end

    function isPairActive = processMortarPair(obj,is,im,elS,elM,k)

      isPairActive = true;

      [xiSlave,xiMaster,detJ] = segmentBasedCouple(obj,is,im,elS,elM);
      if isempty(xiSlave)
        isPairActive = false;
        return
      end

      nTri = numel(detJ);

      for i = 1:nTri
        idTri = (k-1)*obj.maxTriPerPair + i;
        s = MortarSide.slave;
        m = MortarSide.master;
        obj.interfacePairs(idTri,[s m]) = [is im];
        obj.gpCoords{s}(idTri,:,1) = xiSlave(:,1,i);
        obj.gpCoords{s}(idTri,:,2) = xiSlave(:,2,i);
        obj.gpCoords{m}(idTri,:,1) = xiMaster(:,1,i);
        obj.gpCoords{m}(idTri,:,2) = xiMaster(:,2,i);
        obj.detJtri(idTri) = detJ(i);
      end
    end


    function finalizeMortarMaps(obj)
      % remove useless entries after mortar preallocation
      id = obj.detJtri == 0;
      obj.detJtri = obj.detJtri(~id);
      for s = [MortarSide.master,MortarSide.slave]
      obj.gpCoords{s}(id,:,:) = [];
      end
      obj.interfacePairs(id,:) = [];

      obj.numbInterfacePairs = size(obj.interfacePairs,1);
    end

    
    function [xiSlave,xiMaster,dJTri] = segmentBasedCouple(obj,idSlave,idMaster,elemSlave,elemMaster)
      % output: nGx2xnT matrices of reference coordinate in slave and
      % master side to perform segment based integration
      % dJTri: determinant for each pallet
      % nT: numb. of triangles from delaunay triangulation on clipping
      % polygon

      % compute auxiliary plane for integration
      master = MortarSide.master;
      slave = MortarSide.slave;
      elId([slave,master]) = [idSlave,idMaster];
      elem([slave,master]) = [elemSlave,elemMaster];
      xi = cell(2,1);

      ns = zeros(2,1);
      % get number of subsegment ns (in case of higher-order elements)
      for side = [slave,master]
        if elem(side).nNode > 4
          ns(side) = 4;
          elem(side) = elem(side).subQuad;
        else
          ns(side) = 1;
        end
      end


      % initialize output
      nTri = 0;
      [xi{master},xi{slave}] = deal(zeros(numel(obj.wTri),2,obj.maxTriPerPair));
      dJTri = zeros(obj.maxTriPerPair,1);

      % double loop over slave and master subsegments
      for iS = 1:ns(slave)
        if ns(slave) == 1
          [P0,nP] = computeAuxiliaryPlane(obj,elId(slave));
          coordS3D = getElementCoords(elem(slave),elId(slave));
        else
          [P0,nP] = computeAuxiliaryPlane(obj,elId(slave),iS);
          coordS3D =  elem(slave).getSubElementCoords(elId(slave),iS);  
        end

        coordS = pointToSurfaceProjection(P0,nP,coordS3D);

        for iM = 1:ns(master)
          if ns(master) == 1
            coordM3D = getElementCoords(elem(master),elId(master));
          else
            coordM3D = elem(master).getSubElementCoords(elId(master),iM);
          end

          coordM = pointToSurfaceProjection(P0,nP,coordM3D);

          [coordClip,topolClip,isClipValid] = SegmentBasedQuadrature.segmentation(coordS,coordM);

          if ~isClipValid
            continue
          end

          nTriLoc = size(topolClip,1);

          % project gauss points of each triangular cell into slave and
          % master subsegments
          xiSlaveLoc = projectBack(obj,elem(slave),topolClip,coordClip,coordS);
          xiMasterLoc = projectBack(obj,elem(master),topolClip,coordClip,coordM);

          % map subsegment coords to higher order element coords
          if ns(master) > 1
            xiMasterLoc = elem(master).mapsub2ref(xiMasterLoc,iM);
          end
          if ns(slave) > 1
            xiSlaveLoc = elem(slave).mapsub2ref(xiSlaveLoc,iS);
          end

          for iT = 1:nTriLoc
            triVert = coordClip(topolClip(iT,:),:);
            dJTri(nTri+iT) = 2 * Triangle.computeArea(triVert);
          end

          xi{master}(:,:,nTri+1:nTri+nTriLoc) = xiMasterLoc;
          xi{slave}(:,:,nTri+1:nTri+nTriLoc) = xiSlaveLoc;

          nTri = nTri + nTriLoc;
        end
      end

      % cut outputs
      xiMaster = xi{master}(:,:,1:nTri);
      xiSlave = xi{slave}(:,:,1:nTri);
      dJTri = dJTri(1:nTri);

    end
    
    function [P,n] = computeAuxiliaryPlane(obj,el,subID)
      if nargin == 2
        P = obj.grids(MortarSide.slave).surfaces.center(el,:);
        n = obj.grids(MortarSide.slave).surfaces.normal(el,:);
      elseif nargin == 3
        % compute auxiliary plane on subsegment of quad9
        P = elem.computeCentroid(el,subID);
        n = elem.computeNormal(el,elem.centroid,subID);
      end
    end

    function xi = projectBack(obj,elem,topolTri,clipCoord,elemCoord)
      % return reference coordinates in master/slave space for GP in
      % triangle facets after intersection
      % output is a 3D matrices of size nGx2xnTri
      xi = zeros(numel(obj.wTri),2,size(topolTri,1));
      % netwon params
      itMax = 10;
      tol = 1e-9;
      tri = Triangle('gaussOrder',obj.triGaussOrder);                % define reference triangle
      for i = 1:size(topolTri,1)
        coordTri = clipCoord(topolTri(i,:),:);
        coordGPtri = getGPointsLocation(tri,coordTri);
        for g = 1:numel(obj.wTri)
          rhs = (elem.computeBasisF(xi(g,:,i))*elemCoord)' - coordGPtri(g,:)';
          iter = 0;
          while (norm(rhs,2) > tol) && (iter < itMax)
            iter = iter+1;
            J = (elem.computeDerBasisF(xi(g,:,i))*elemCoord)';
            dxi = J\(-rhs);
            xi(g,:,i) = xi(g,:,i) + dxi';
            rhs = (elem.computeBasisF(xi(g,:,i))*elemCoord)' - coordGPtri(g,:)';
          end
          fl = elem.checkInRange(xi(g,:,i));
          assert(fl,['GP %i in triangle %i out of range of coordinates of' ...
            'reference space']);
        end
      end
    end


    function dJw = getIntegrationWeights(obj,kPair)

      dJ = obj.detJtri(kPair);
      dJw = dJ*obj.wTri;

    end

    function gpCoords = getSlaveGPCoords(obj,kPair)
      gpCoords = obj.gpCoords{MortarSide.slave}(kPair,:,:);
      gpCoords = reshape(gpCoords,[],2);
    end

    function gpCoords = getMasterGPCoords(obj,kPair)
      gpCoords = obj.gpCoords{MortarSide.master}(kPair,:,:);
      gpCoords = reshape(gpCoords,[],2);
    end


  end

  methods (Static)


    function isValid = isClipValid(clip)
      tol = 1e-4;
      clip = round(clip/tol)*tol;
      clipUnique = unique(clip,'rows');
      isValid = true;
      if size(clipUnique,1) < 3
        isValid = false;
      end

    end

    function [polyClip,topolClip,isClipValid] = segmentation(poly1,poly2)

      polyClip = clipPolygon(poly1,poly2);

      topolClip = [];

      if isempty(polyClip)
        isClipValid = false;
        return
      end

      %fan triangulation
      nV = size(polyClip,1);
      isClipValid = true;
      topolClip = [ ones(nV-2,1), (2:nV-1)', (3:nV)' ];

    end

  
  end
end

