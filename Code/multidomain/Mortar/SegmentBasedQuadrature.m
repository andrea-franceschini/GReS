classdef SegmentBasedQuadrature < handle
  % Implement utilities to perform segment based integration given a pair
  % of master and slave element in contact
  % Can be equipped with the RBF scheme to cheaply sort out elements that
  % are not really connected

  % REFS: Puso,2004, A mortar segment-to-segment contact method for large
  % deformation solid mechanics
  
  properties
    msh       % instance of InterfaceMesh class
    mortar
    nGtri
    elems     % provisional elems instances for each master/slave pair
    dJw
    nTri
    maxTri = 30;
  end
  
  methods
    function obj = SegmentBasedQuadrature(mortar,nGtri)
      obj.mortar = mortar;
      obj.msh = obj.mortar.mesh.msh;
      obj.nGtri = nGtri;
    end

    function [Ns,Nm,Nmult,NbubSlave,NbubMaster] = getMortarBasisFunctions(obj,is,im)
      assert(nargout >=3, 'Not enough output arguments');
      [xiSlave,xiMaster,obj.dJw] = segmentBasedCouple(obj,is,im);
      if isempty(xiSlave)
        [Ns,Nm,Nmult,NbubSlave,NbubMaster] = deal([]);
        return
      end
      elemSlave = obj.mortar.getElem(2,is);
      elemMaster = obj.mortar.getElem(1,im);

      Nm = zeros(obj.nGtri, elemMaster.nNode,obj.nTri);
      Ns = zeros(obj.nGtri,elemSlave.nNode,obj.nTri);
      if strcmp(obj.mortar.multiplierType,'P0')
        d2 = 1;
      else
        d2 = elemSlave.nNode;
      end
      Nmult = zeros(obj.nGtri,d2,obj.nTri);

      if nargout > 3
        [NbubSlave,NbubMaster] = deal(zeros(obj.nGtri,1,obj.nTri));
      end

      for i = 1:obj.nTri
        xiM = xiMaster(:,:,i);
        xiS = xiSlave(:,:,i);
        Nm(:,:,i) = elemMaster.computeBasisF(xiM);
        Ns(:,:,i) = elemSlave.computeBasisF(xiS);
        Nmult(:,:,i) = computeMultiplierBasisF(obj.mortar,is,Ns(:,:,i));
        if nargout > 3
          NbubSlave(:,:,i) = elemSlave.computeBubbleBasisF(xiS);
          NbubMaster(:,:,i) = elemMaster.computeBubbleBasisF(xiM);
        end
      end
    end


    function mat = integrate(obj,func,varargin)
      % perform segment based integration
      % func: handle to the integrand
      % varargin: list of integrand inputs
      % each integrand is a 4D matrix of size: ndofRow,ndofcol,nG,nTri

      % check input
      assert(nargin(func)==numel(varargin),['Number of specified input (%i)' ...
        'not matching the integrand input (%i)'],numel(varargin),nargin(func));
      size4 = cellfun(@(x) size(x, 4), varargin);
      if ~all(size4 == obj.nTri)
        error('All inputs must have the same size along dimension 3.');
      end

      % segment based loop over triangular pallets

      for i = 1:obj.nTri
        args = cellfun(@(x) x(:, :, :, i), varargin, 'UniformOutput', false);
        matloc = func(args{:});
        matloc = matloc.*reshape(obj.dJw(:,i),1,1,[]);
        matloc = sum(matloc,3);
        if i ==1
          mat = matloc;
        else
          mat = mat + matloc;
        end
      end
    end 
    
    function [xiSlave,xiMaster,dJwTri] = segmentBasedCouple(obj,elSlave,elMaster)
      % output: nGx2xnT matrices of reference coordinate in slave and
      % master side to perform segment based integration
      % dJwTri: weighed jacobian determinant for each pallet. size nGxnTri
      % nT: numb. of triangles from delaunay triangulation on clipping
      % polygon

      % compute auxiliary plane for integration
      obj.elems = {getElem(obj.mortar,1,elMaster),...
        getElem(obj.mortar,2,elSlave)};

      % get number of susegment (in case of higher-order elements)
      if obj.elems{1}.nNode > 4
        ns1 = 4;
        elemS = createSubElement(obj.elems{2});
      else
        ns1 = 1;
        elemS = obj.elems{2};
      end

      if obj.elems{2}.nNode > 4
        ns2 = 4;
        elemM = createSubElement(obj.elems{1});
      else
        ns2 = 1;
        elemM = obj.elems{1};
      end

      % initialize output
      obj.nTri = 0;
      [xiMaster,xiSlave] = deal(zeros(obj.nGtri,2,obj.maxTri));
      dJwTri = zeros(obj.nGtri,obj.maxTri);

      % double loop over slave and master subsegments
      for iS = 1:ns2
        if ns2 == 1
          [P0,nP] = computeAuxiliaryPlane(obj,elSlave);
          coordS3D = FEM.getElementCoords(obj.elems{2},elSlave);
        else
          [P0,nP] = computeAuxiliaryPlane(obj,elSlave,iS);
          coordS3D =  obj.elems{2}.getSubElementCoords(elSlave,iS);  
        end

        coordS = projectNodes(obj,P0,nP,coordS3D);

        for iM = 1:ns1
          if ns1==1
            coordM3D = FEM.getElementCoords(obj.elems{1},elMaster);
          else
            coordM3D = obj.elems{1}.getSubElementCoords(elMaster,iM);
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
          topolClip = delaunay(coordClip(:,1),coordClip(:,2));
          nTriLoc = size(topolClip,1);

          % project gauss points of each triangular cell into slave and
          % master susegments
          xiSlaveLoc = projectBack(obj,elemS,topolClip,coordClip,coordS);
          xiMasterLoc = projectBack(obj,elemM,topolClip,coordClip,coordM);

          % map subsegment coords to higher order element coords
          if ns1 > 1
            xiMasterLoc = obj.elems{1}.mapsub2ref(xiMasterLoc,iM);
          end
          if ns2 > 1
            xiSlaveLoc = obj.elems{1}.mapsub2ref(xiSlaveLoc,iS);
          end

          % compute determinant in triangles
          tri = Triangle(obj.nGtri);
          for iT = 1:nTriLoc
            triVert = coordClip(topolClip(iT,:),:);
            dJwTri(:,obj.nTri+iT) = getDerBasisFAndDet(tri,triVert);
          end

          xiMaster(:,:,obj.nTri+1:obj.nTri+nTriLoc) = xiMasterLoc;
          xiSlave(:,:,obj.nTri+1:obj.nTri+nTriLoc) = xiSlaveLoc;

          obj.nTri = obj.nTri + nTriLoc;
        end
      end

      % cut outputs
      if obj.nTri == 0
        [xiSlave,xiMaster,dJwTri] = deal([]);
      end
      xiMaster = xiMaster(:,:,1:obj.nTri);
      xiSlave = xiSlave(:,:,1:obj.nTri);
      dJwTri = dJwTri(:,1:obj.nTri);

    end
    
    function [P,n] = computeAuxiliaryPlane(obj,el,subID)
      if nargin == 2
        P = obj.msh(2).surfaceCentroid(el,:);
        n = obj.elems{2}.computeNormal(el,obj.elems{2}.centroid);
      elseif nargin == 3
        % compute auxiliary plane on subsegment of quad9
        P = obj.elems{2}.computeCentroid(el,subID);
        n = obj.elems{2}.computeNormal(el,obj.elems{2}.centroid,subID);
      end
    end

    function xi = projectBack(obj,elem,topolTri,clipCoord,elemCoord)
      % return reference coordinates in master/slave space for GP in
      % triangle facets after intersection
      % output is a 3D matrices of size nGx2xnTri
      xi = zeros(obj.nGtri,2,size(topolTri,1));
      % netwon params
      itMax = 10;
      tol = 1e-9;
      tri = Triangle(obj.nGtri);                % define reference triangle
      for i = 1:size(topolTri,1)
        coordTri = clipCoord(topolTri(i,:),:);
        coordGPtri = getGPointsLocation(tri,coordTri);
        for g = 1:obj.nGtri
          rhs = (elem.computeBasisF(xi(g,:,i))*elemCoord)' - coordGPtri(g,:)';
          iter = 0;
          while (norm(rhs,2) > tol) && (iter < itMax)
            iter = iter+1;
            J = (elem.computeDerBasisF(xi(g,:,i))*elemCoord)';
            dxi = J\(-rhs);
            xi(g,:,i) = xi(g,:,i) + dxi';
            rhs = (elem.computeBasisF(xi(g,:,i))*elemCoord)' - coordGPtri(g,:)';
          end
          fl = FEM.checkInRange(elem,xi(g,:,i));
          %assert(fl,['Error for GP %i in triangle %i out of range of ' ...
          %  'reference space coordinates.'],g,i);
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
  end
end

