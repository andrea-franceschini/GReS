classdef SegmentBasedQuadrature < MortarQuadrature
  % Implement utilities to perform segment based integration given a pair
  % of master and slave element in contact
  % Can be equipped with the RBF scheme to cheaply sort out elements that
  % are not really connected

  % REFS: Puso,2004, A mortar segment-to-segment contact method for large
  % deformation solid mechanics
  
  properties
    detJtri
    ngTri
  end

  properties (Access = private)
    elems     % provisional elems instances for each master/slave pair
    wTri      % store weight for gauss integration on triangles
    maxTriPerPair
  end
  
  methods
    function obj = SegmentBasedQuadrature(interface,multType,input)

      obj@MortarQuadrature(interface,multType,input);
      obj.ngTri = getXMLData(input,[],"nGP");

      gaussTri = Gauss(5,obj.ngTri); % 5 is the vtk type of triangles
      obj.wTri = gaussTri.weight;

      isQuadratic = [~isempty(obj.elements(1).getElement(28));
                     ~isempty(obj.elements(2).getElement(28))];

      if sum(isQuadratic)==0 % no quadratic elements
        obj.maxTriPerPair = 6;
      elseif sum(isQuadratic)==1
        obj.maxTriPerPair = 24;
      elseif sum(isQuadratic)==2
        obj.maxTriPerPair = 96;
      end


    end

    function processMortarPairs(obj)

      % initialize the maps to store mortar quadrature infos
      nConnections = nnz(obj.interface.interfMesh.elemConnectivity);
      totTri = nConnections*obj.maxTriPerPair;
      obj.gpCoords = {zeros(totTri,obj.ngTri,2);
        zeros(totTri,obj.ngTri,2)};
      obj.interfacePairs = zeros(totTri,2);
      obj.detJtri = zeros(totTri,1);

      nM = full(sum(obj.interface.interfMesh.elemConnectivity,1));
      nM = [0 cumsum(nM)];

      for is = 1:obj.msh(2).nSurfaces

        imList = find(obj.interface.interfMesh.elemConnectivity(:,is));

        for j = 1:numel(imList)
          im = imList(j);
          k = nM(is)+ j;
          % k (global index to write without race conditions)
          processMortarPair(obj,is,im,k);
        end
      end

      finalizeMortarMaps(obj);
      computeAreaSlave(obj);

    end

    function isPairActive = processMortarPair(obj,is,im,k)

      isPairActive = true;

      [xiSlave,xiMaster,detJ] = segmentBasedCouple(obj,is,im);
      if isempty(xiSlave)
        isPairActive = false;
        return
      end

      nTri = numel(detJ);

      for i = 1:nTri
        idTri = (k-1)*obj.maxTriPerPair + i;
        obj.interfacePairs(idTri,:) = [is im];
        obj.gpCoords{2}(idTri,:,1) = xiSlave(:,1,i);
        obj.gpCoords{2}(idTri,:,2) = xiSlave(:,2,i);
        obj.gpCoords{1}(idTri,:,1) = xiMaster(:,1,i);
        obj.gpCoords{1}(idTri,:,2) = xiMaster(:,2,i);
        obj.detJtri(idTri) = detJ(i);
      end
    end


    function finalizeMortarMaps(obj)
      % remove useless entries after mortar preallocation
      id = obj.detJtri == 0;
      obj.detJtri = obj.detJtri(~id);
      for i = 1:2
      obj.gpCoords{i}(id,:,:) = [];
      end
      obj.interfacePairs(id,:) = [];

      obj.numbInterfacePairs = size(obj.interfacePairs,1);
    end

    
    function [xiSlave,xiMaster,dJTri] = segmentBasedCouple(obj,elSlave,elMaster)
      % output: nGx2xnT matrices of reference coordinate in slave and
      % master side to perform segment based integration
      % dJTri: determinant for each pallet
      % nT: numb. of triangles from delaunay triangulation on clipping
      % polygon

      % compute auxiliary plane for integration
      obj.elems = [getElem(obj,1,elMaster),...
                   getElem(obj,2,elSlave)];

      % get number of susegment (in case of higher-order elements)
      if obj.elems(1).nNode > 4
        ns1 = 4;
        elemS = createSubElement(obj.elems(1));
      else
        ns1 = 1;
        elemS = obj.elems(2);
      end

      if obj.elems(2).nNode > 4
        ns2 = 4;
        elemM = createSubElement(obj.elems(2));
      else
        ns2 = 1;
        elemM = obj.elems(2);
      end

      % initialize output
      nTri = 0;
      [xiMaster,xiSlave] = deal(zeros(obj.ngTri,2,obj.maxTriPerPair));
      dJTri = zeros(obj.maxTriPerPair,1);

      % double loop over slave and master subsegments
      for iS = 1:ns2
        if ns2 == 1
          [P0,nP] = computeAuxiliaryPlane(obj,elSlave);
          coordS3D = FEM.getElementCoords(obj.elems(2),elSlave);
        else
          [P0,nP] = computeAuxiliaryPlane(obj,elSlave,iS);
          coordS3D =  obj.elems(2).getSubElementCoords(elSlave,iS);  
        end

        coordS = projectNodes(obj,P0,nP,coordS3D);

        for iM = 1:ns1
          if ns1==1
            coordM3D = FEM.getElementCoords(obj.elems(1),elMaster);
          else
            coordM3D = obj.elems(1).getSubElementCoords(elMaster,iM);
          end

          coordM = projectNodes(obj,P0,nP,coordM3D);

          % compute 2D coordinates of projected nodes on the aux. plane
          % obtain intersection of polygon
          [clipX,clipY] = polyclip(coordS(:,1),coordS(:,2),coordM(:,1),coordM(:,2),1);
          if numel(clipX)==0
            % no intersection
            continue
          end
          assert(numel(clipX)==1,['Non unique clip polygon for master/slave pair' ...
            ' %i/%i \n'],elMaster,elSlave)


          % perform delaunay triangulation on clip polygon
          % assumption: only one clip polygon results from intersection
          coordClip = [clipX{:} clipY{:}];

          if ~SegmentBasedQuadrature.isClipValid(coordClip)
            % skip if the polygon is degenerate or has very small area
            continue
          end
          topolClip = delaunay(coordClip(:,1),coordClip(:,2));
          nTriLoc = size(topolClip,1);

          % project gauss points of each triangular cell into slave and
          % master susegments
          xiSlaveLoc = projectBack(obj,elemS,topolClip,coordClip,coordS);
          xiMasterLoc = projectBack(obj,elemM,topolClip,coordClip,coordM);

          % map subsegment coords to higher order element coords
          if ns1 > 1
            xiMasterLoc = obj.elems(1).mapsub2ref(xiMasterLoc,iM);
          end
          if ns2 > 1
            xiSlaveLoc = obj.elems(1).mapsub2ref(xiSlaveLoc,iS);
          end

          for iT = 1:nTriLoc
            triVert = coordClip(topolClip(iT,:),:);
            dJTri(nTri+iT) = 2 * Triangle.computeArea(triVert);
          end

          xiMaster(:,:,nTri+1:nTri+nTriLoc) = xiMasterLoc;
          xiSlave(:,:,nTri+1:nTri+nTriLoc) = xiSlaveLoc;

          nTri = nTri + nTriLoc;
        end
      end

      % cut outputs
      if nTri == 0
        [xiSlave,xiMaster,dJTri] = deal([]);
      end
      xiMaster = xiMaster(:,:,1:nTri);
      xiSlave = xiSlave(:,:,1:nTri);
      dJTri = dJTri(1:nTri);

    end
    
    function [P,n] = computeAuxiliaryPlane(obj,el,subID)
      if nargin == 2
        P = obj.msh(2).surfaceCentroid(el,:);
        n = obj.elems(2).computeNormal(el,obj.elems(2).centroid);
      elseif nargin == 3
        % compute auxiliary plane on subsegment of quad9
        P = obj.elems(2).computeCentroid(el,subID);
        n = obj.elems(2).computeNormal(el,obj.elems(2).centroid,subID);
      end
    end

    function xi = projectBack(obj,elem,topolTri,clipCoord,elemCoord)
      % return reference coordinates in master/slave space for GP in
      % triangle facets after intersection
      % output is a 3D matrices of size nGx2xnTri
      xi = zeros(obj.ngTri,2,size(topolTri,1));
      % netwon params
      itMax = 10;
      tol = 1e-9;
      tri = Triangle(obj.ngTri);                % define reference triangle
      for i = 1:size(topolTri,1)
        coordTri = clipCoord(topolTri(i,:),:);
        coordGPtri = getGPointsLocation(tri,coordTri);
        for g = 1:obj.ngTri
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

    function x = projectNodes(obj,P,n,coord)
      % Project nodes of triangle pair into auxiliary plane
      % get 2D direction of plane
      % Choose arbitrary vector not parallel to n
      if abs(n(1)) < 0.9
        temp = [1; 0; 0];
      else
        temp = [0; 1; 0];
      end
      % First direction in the plane
      d1 = cross(n, temp);
      d1 = d1 / norm(d1);
      % Second direction in the plane
      d2 = cross(n, d1);

      c = coord;
      cn = (c - P)*n;
      nN = size(c,1);
      projC = c - repmat(n',nN,1).*cn;
      % project in 2D
      x = (projC - P)*[d1 d2];

    end


    function dJw = getIntegrationWeights(obj,kPair)

      dJ = obj.detJtri(kPair);
      dJw = dJ*obj.wTri;

    end

    function gpCoords = getSlaveGPCoords(obj,kPair)
      gpCoords = obj.gpCoords{2}(kPair,:,:);
      gpCoords = reshape(gpCoords,[],2);
    end

    function gpCoords = getMasterGPCoords(obj,kPair)
      gpCoords = obj.gpCoords{1}(kPair,:,:);
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

  end
end

