classdef Mortar2D < handle
    % Class implementing algorithms for mortar method with 1D interfaces 
    properties
        mshMaster
        mshSlave
        masterTopol
        slaveTopol
        masterCoord
        slaveCoord
        nElMaster
        nElSlave
        nNodesMaster
        nNodesSlave
        elemConnectivity
        degree 
        elemsMaster
        elemsSlave
        nodesMaster
        nodesSlave
        Dmat
        nSMat
        nMMat
        nNelem
    end
    
    methods
        function obj = Mortar2D(deg,varargin)
            % Constructor for the class
            % INPUT:
            % empty varargin (later manual assignment of
            % coordinates and topology of master and slave
            % non empty varargin: meshMaster, id master interface,
            % meshSlave, id slave interface
            obj.degree = deg;
            if deg == 1
                obj.nNelem = 2;
            else
                obj.nNelem = 3;
            end
            if  nargin == 5
                obj.mshMaster = varargin{1};
                obj.mshSlave = varargin{3};
                masterTop = varargin{1}.edges;
                obj.masterCoord = varargin{1}.coordinates;
                obj.masterTopol = masterTop(varargin{1}.edgeTag == varargin{2}, :);
                slaveTop = varargin{3}.edges;
                obj.slaveCoord = varargin{3}.coordinates;
                obj.slaveTopol = slaveTop(varargin{3}.edgeTag == varargin{4}, :);
            elseif nargin == 6 && strcmp(varargin{1},'set')
                obj.masterTopol = varargin{2};
                obj.slaveTopol = varargin{3};
                obj.masterCoord = varargin{4};
                obj.slaveCoord = varargin{5};
            else
                error('Wrong number of inputs for class Mortar1D')
            end
            % obj.nodesMaster = unique(obj.masterTopol);
            % obj.nodesSlave = unique(obj.slaveTopol);
            obj.nElMaster = size(obj.masterTopol,1);
            obj.nElSlave = size(obj.slaveTopol,1);
            obj.nNodesMaster = size(obj.masterCoord,1);
            obj.nNodesSlave = size(obj.slaveCoord,1);
            obj.elemConnectivity = computeElementConnectivity(obj);
            % connected master and slave surfaces (find non empty rows and
            % columns of connectivity matrix)
            idM = sum(obj.elemConnectivity,2) > 0;
            idS = sum(obj.elemConnectivity,1) > 0;
            % get list of nodes that actually belong to elements in contact
            obj.elemsMaster = find(idM);
            obj.elemsSlave = find(idS);
            obj.nodesMaster = unique(obj.masterTopol(idM,:));
            obj.nodesSlave = unique(obj.slaveTopol(idS,:));
            getMatricesSize(obj);
            %computeSlaveMatrix(obj);
        end

        function obj = computeSlaveMatrix(obj)
            D = zeros(obj.nSMat,obj.nSMat);
            if obj.degree < 2
                for sID = 1:obj.nElSlave
                    % elem length
                    i1 = obj.slaveTopol(sID,1);
                    i2 = obj.slaveTopol(sID,2);
                    h = norm(obj.slaveCoord(i1,:)-obj.slaveCoord(i2,:));
                    Dloc = (h/6)*[2 1; 1 2];
                    D([i1 i2], [i1 i2]) = D([i1 i2], [i1 i2]) + Dloc;
                end
            else
                g = Gauss(12,2,1);
                gpRef = g.coord;
                gpWeight = g.weight;
                for sID = 1:size(obj.slaveTopol,1)
                    % elem length
                    i1 = obj.slaveTopol(sID,1);
                    i2 = obj.slaveTopol(sID,2);
                    i3 = obj.slaveTopol(sID,3);
                    h = norm(obj.slaveCoord(i1,:)-obj.slaveCoord(i3,:));
                    N = computeQuadraticSF(gpRef);
                    Dloc = 0.5*h*N'*(N.*gpWeight);
                    D([i1 i2 i3], [i1 i2 i3]) = D([i1 i2 i3], [i1 i2 i3]) + Dloc;
                end
            end
            % extract only mass matrix entries belonging to the interface
            obj.Dmat = D(obj.nodesSlave, obj.nodesSlave);
        end



        function [D,M,varargout] = computeMortarSegmentBased(obj,nGP,mult_type)
            if length(unique(obj.slaveCoord(obj.nodesSlave,2))) ~= 1
                if nargout > 2
                    varargout{1} = 0;
                end
                return
            end
            M = zeros(obj.nSMat,obj.nMMat);
            D = zeros(obj.nSMat,obj.nSMat);
            for i = 1:obj.nElMaster
                % get element of the support
                a = obj.masterCoord(obj.masterTopol(i,1),:);
                b = obj.masterCoord(obj.masterTopol(i,2),:);
                idMaster = obj.masterTopol(i,1:end);
                % get shape function values and interpolation points coordinate
                slave_elems = find(obj.elemConnectivity(i,:));
                % loop trough connected slave elements
                for j = slave_elems
                  % get nodes
                  g = Gauss(obj.mshSlave.edgeVTKType(j),nGP);
                  gpRef = g.coord;
                  gpWeight = g.weight;
                  c = obj.slaveCoord(obj.slaveTopol(j,1),:);
                  d = obj.slaveCoord(obj.slaveTopol(j,2),:);
                  % find intersection between master and slave (just looking at the
                  % x-projection)
                    tmp = sort([a(1) b(1) c(1) d(1)]);
                    i1 = [tmp(2) 0];
                    i2 = [tmp(3) 0];
                    % get GP in the intersection element
                    gPts = ref2nod(gpRef, i1, i2);
                    % compute basis function on master element points
                    gMaster = nod2ref(gPts, a, b);
                    Nmaster = computeBasis1D(gMaster,obj.degree);
                    % compute basis function on slave elements points
                    gSlave = nod2ref(gPts, c, d);
                    Nslave = computeBasis1D(gSlave,obj.degree);
                    switch mult_type
                      case 'standard'
                        Nmult = Nslave;
                      case 'dual'
                        % compute dual basis function
                        Nmult = computeDualBasisF(obj,gSlave);
                    end
                    Mloc = Nmult'*(Nmaster.*gpWeight);
                    Dloc = Nmult'*(Nslave.*gpWeight);
                    Mloc = Mloc*(0.5*(i2(1)-i1(1)));
                    Dloc = Dloc*(0.5*(i2(1)-i1(1)));
                    idSlave = obj.slaveTopol(j,1:end);
                    M(idSlave, idMaster) = M(idSlave, idMaster) + Mloc;
                    D(idSlave,idSlave) = D(idSlave,idSlave) + Dloc;
                end
            end
            M = M(obj.nodesSlave, obj.nodesMaster);
            D = D(obj.nodesSlave, obj.nodesSlave);
            if nargout == 3
                varargout{1} = D\M;
            end
        end

        function [D,M,varargout] = computeMortarElementBased(obj,nGP,mult_type)
            n = computeNodalNormal(obj);
            M = zeros(obj.nSMat,obj.nMMat);
            D = zeros(obj.nSMat,obj.nSMat);
            nN = 2; % numb. nodes per element
            if obj.degree > 1
                nN = 3;
            end
            for i = 1:obj.nElSlave
                % get nodes
                g = Gauss(obj.mshSlave.edgeVTKType(i),nGP);
                gpRef = g.coord;
                gpWeight = g.weight;
                gpPos = gpRef;
                gpW = gpWeight;
                s1 = obj.slaveCoord(obj.slaveTopol(i,1),:);
                s2 = obj.slaveCoord(obj.slaveTopol(i,2),:);
                idSlave = obj.slaveTopol(i,1:nN);
                xSlave = ref2nod(gpPos, s1, s2);
                master_elems = find(obj.elemConnectivity(:,i));
                %id = true(length(gpRef),1);
                for m = master_elems'
                   [xiM] = obj.projectGP(m,i,gpPos,n,xSlave);
                   id = all([xiM >= -1, xiM <= 1],2);
                   if any(id)
                      NMaster = computeBasis1D(xiM(id),obj.degree);
                      NSlave = computeBasis1D(gpPos(id),obj.degree);
                      switch mult_type
                        case 'standard'
                          Nmult = NSlave;
                        case 'dual'
                          % compute dual basis function
                          Nmult = computeDualBasisF(obj,gpPos(id));
                      end
                      Mloc = Nmult'*(NMaster.*gpW(id));
                      Mloc = Mloc*(0.5*sqrt((s1(1)-s2(1))^2 + (s1(2)-s2(2))^2));
                      Dloc = Nmult'*(NSlave.*gpW(id));
                      Dloc = Dloc*(0.5*sqrt((s1(1)-s2(1))^2 + (s1(2)-s2(2))^2));
                      idMaster = obj.masterTopol(m,:);
                      M(idSlave, idMaster) = M(idSlave, idMaster) + Mloc;
                      D(idSlave, idSlave) = D(idSlave, idSlave) + Dloc;
                      % sort out Points already projected
                      gpPos = gpPos(~id);
                      gpW = gpW(~id);
                      xSlave = xSlave(~id,:);
                   end
                end
            end
            D = D(obj.nodesSlave,obj.nodesSlave);
            M = M(obj.nodesSlave, obj.nodesMaster);
            if nargout == 3
               varargout{1} = D\M;
            end
        end

        function [D,M,varargout] = computeMortarRBF(obj,nGP,nInt,type,mult_type)
           % type: family of Radial Basis function to use
           % mult_type:
           tol = 1e-3;
           M = zeros(obj.nSMat,obj.nMMat);
           D = zeros(obj.nSMat,obj.nSMat);
           [wFMat,w1Mat,ptsIntMat] = getWeights(obj,nInt,type);
           % Loop trough slave elements
           for j = 1:obj.nElSlave
             if j>1 && ~all(id)
               %error('%i GP not projected for element %i',sum(~id),j-1)
             end
             g = Gauss(obj.mshSlave.edgeVTKType(j),nGP);
             gpRef = g.coord;
             gpW = g.weight;
             % nodes of slave element
             s1 = obj.slaveCoord(obj.slaveTopol(j,1),:);
             s2 = obj.slaveCoord(obj.slaveTopol(j,2),:);
              h = sqrt((s1(1)-s2(1))^2 + (s1(2)-s2(2))^2);
              idSlave = obj.slaveTopol(j,1:obj.nNelem);
              % Real position of Gauss points
              ptsGauss = ref2nod(gpRef, s1, s2);
              % Get connected master elements
              master_elems = find(obj.elemConnectivity(:,j));
              NSlave = computeBasis1D(gpRef,obj.degree);
              proj_GP = 1:g.nNode; % list of GP to project
              switch mult_type
                 case 'standard'
                    NSlaveMult = NSlave;
                 case 'dual'
                    % compute dual basis function
                    NSlaveMult = computeDualBasisF(obj,gpRef);
              end
              % Loop trough master elements and store projected master
              % basis functions
              for jm = master_elems'
                 ptsInt = ptsIntMat(:,repNum(2,jm));
                 idMaster = obj.masterTopol(jm,1:obj.nNelem);
                 fiNM = obj.computeRBFfiNM(ptsInt,ptsGauss,type);
                 %NMaster = fiNM*wFMat(:,getWeightsID(obj,jm));
                 switch obj.degree
                    case 1
                       NMaster = (fiNM*wFMat(:,repNum(obj.nNelem,jm)))./(fiNM*w1Mat(:,jm));
                       NSupp = NMaster(:,1);
                    case 2
                       Ntmp = (fiNM*wFMat(:,repNum(obj.nNelem+1,jm)))./(fiNM*w1Mat(:,jm));
                       NMaster = Ntmp(:,1:obj.nNelem);
                       NSupp = Ntmp(:,end);
                 end
                 id = all([NSupp >= 0 - tol, NSupp <= 1 + tol],2);
                 proj_GP = proj_GP(~id);
                 if any(id)
                    NMaster = NMaster(id,:);
                    Mloc = NSlaveMult(id,:)'*(NMaster.*gpW(id));
                    Mloc = Mloc*(0.5*h);
                    Dloc = NSlaveMult(id,:)'*(NSlave(id,:).*gpW(id));
                    Dloc = Dloc*(0.5*h);
                    M(idSlave, idMaster) = M(idSlave, idMaster) + Mloc;
                    D(idSlave, idSlave) = D(idSlave, idSlave) + Dloc;
                    % sort out Points already projected
                    gpRef = gpRef(~id);
                    gpW = gpW(~id);
                    ptsGauss = ptsGauss(~id,:);
                    NSlave = NSlave(~id,:);
                    NSlaveMult = NSlaveMult(~id,:);
                 end
              end
           end
              M = M(obj.nodesSlave, obj.nodesMaster);
              D = D(obj.nodesSlave, obj.nodesSlave);
              D(abs(D)<eps) = 0;
           %end
           %t = toc;
           if nargout == 3
              varargout{1} = t;
           elseif nargout == 4
              varargout{1} = t;
              varargout{2} = D\M;
           end
        end


        function [D,M,varargout] = computeMortarRBF_new(obj,nGP,nInt,type,mult_type)
           % type: family of Radial Basis function to use
           % mult_type:
           tic
           g = Gauss(12,nGP,1);
           M = zeros(obj.nSMat,obj.nMMat);
           D = zeros(obj.nSMat,obj.nSMat);
           [wFMat,w1Mat,ptsIntMat] = getWeights(obj,nInt,type);
           % Loop trough slave elements
           for j = 1:obj.nElSlave
              Mtmp = M;
              Dtmp = D;
              gpPos = g.coord;
              gpW = g.weight;
              % nodes of slave element
              s1 = obj.slaveCoord(obj.slaveTopol(j,1),:);
              s2 = obj.slaveCoord(obj.slaveTopol(j,end),:);
              h = sqrt((s1(1)-s2(1))^2 + (s1(2)-s2(2))^2);
              idSlave = obj.slaveTopol(j,1:obj.nNelem);
              % Real position of Gauss points
              ptsGauss = ref2nod(gpPos, s1, s2);
              % Get connected master elements
              master_elems = find(obj.elemConnectivity(:,j));
              NSlave = computeBasis1D(gpPos,obj.degree);
              proj_GP = 1:g.nNode; % list of GP to project
              switch mult_type
                 case 'standard'
                    NSlaveMult = NSlave;
                 case 'dual'
                    % compute dual basis function
                    NSlaveMult = computeDualBasisF(obj,gpW,gpPos,h);
              end
              % Loop trough master elements
              % store projected master basis functions
              % find projection pattern and possibly non projecting slave
              % GP (boundary elements situation)
              m=0;
              NMaster = zeros(nGP,2);
              for jm = master_elems'
                 m = m+1;
                 proj_elem = false(nueml(master_elems),nGP);
                 ptsInt = ptsIntMat(:,repNum(2,jm));
                 idMaster = obj.masterTopol(jm,1:obj.nNelem);
                 fiNM = obj.computeRBFfiNM(ptsInt,ptsGauss,type);
                 %NMaster = fiNM*wFMat(:,getWeightsID(obj,jm));
                 switch obj.degree
                    case 1
                       NM = (fiNM*wFMat(:,repNum(obj.nNelem,jm)))./(fiNM*w1Mat(:,jm));
                       NSupp = NM(:,1);
                    case 2
                       Ntmp = (fiNM*wFMat(:,repNum(obj.nNelem+1,jm)))./(fiNM*w1Mat(:,jm));
                       NM = Ntmp(:,1:obj.nNelem);
                       NSupp = Ntmp(:,end);
                 end
                 id = all([NSupp >= 0, NSupp <= 1],2);
                 NMaster(id,:) = NM(id);
                 proj_GP = proj_GP(~id);
                 proj_elem(m,proj_GP(id)) = 1;   
              end
                 if any(id)
                    NMaster = NMaster(id,:);
                    Mloc = NSlaveMult(id,:)'*(NMaster.*gpW(id));
                    Mloc = Mloc*(0.5*h);
                    Dloc = NSlaveMult(id,:)'*(NSlave(id,:).*gpW(id));
                    Dloc = Dloc*(0.5*h);
                    M(idSlave, idMaster) = M(idSlave, idMaster) + Mloc;
                    D(idSlave, idSlave) = D(idSlave, idSlave) + Dloc;
                    % sort out Points already projected
                    gpPos = gpPos(~id);
                    gpW = gpW(~id);
                    ptsGauss = ptsGauss(~id,:);
                    NSlave = NSlave(~id,:);
                    NSlaveMult = NSlaveMult(~id,:);
                 end
              end
              if ~isempty(proj_GP)
                 fprintf('BOUNDARY SLAVE DETECTED! \n')
                 M = Mtmp;
                 D = Dtmp;
                 gpPos = g.coord;
                 gpW = g.weight;
                 % remove gp not projected
                 gpPos = gpPos(~ismember(1:nGP,proj_GP));
                 gpW = gpW(~ismember(1:nGP,proj_GP));
                 % Real position of Gauss points
                 ptsGauss = ref2nod(gpPos, s1, s2);
                 NSlave = computeBasis1D(gpPos,obj.degree);
                 switch mult_type
                    case 'standard'
                       NSlaveMult = NSlave;
                    case 'dual'
                       % recompute dual basis function only on active GP
                       NSlaveMult = computeDualBasisF(obj,gpW,gpPos,h);
                 end
                 % Loop troug master elements
                 for jm = master_elems'
                    ptsInt = ptsIntMat(:,repNum(2,jm));
                    idMaster = obj.masterTopol(jm,1:obj.nNelem);
                    fiNM = obj.computeRBFfiNM(ptsInt,ptsGauss,type);
                    %NMaster = fiNM*wFMat(:,getWeightsID(obj,jm));
                    switch obj.degree
                       case 1
                          NMaster = (fiNM*wFMat(:,repNum(obj.nNelem,jm)))./(fiNM*w1Mat(:,jm));
                          NSupp = NMaster(:,1);
                       case 2
                          Ntmp = (fiNM*wFMat(:,repNum(obj.nNelem+1,jm)))./(fiNM*w1Mat(:,jm));
                          NMaster = Ntmp(:,1:obj.nNelem);
                          NSupp = Ntmp(:,end);
                    end
                    id = all([NSupp >= 0, NSupp <= 1],2);
                    proj_GP = proj_GP(~id);
                    if any(id)
                       NMaster = NMaster(id,:);
                       Mloc = NSlaveMult(id,:)'*(NMaster.*gpW(id));
                       Mloc = Mloc*(0.5*h);
                       Dloc = NSlaveMult(id,:)'*(NSlave(id,:).*gpW(id));
                       Dloc = Dloc*(0.5*h);
                       M(idSlave, idMaster) = M(idSlave, idMaster) + Mloc;
                       D(idSlave, idSlave) = D(idSlave, idSlave) + Dloc;
                       % sort out Points already projected
                       gpPos = gpPos(~id);
                       gpW = gpW(~id);
                       ptsGauss = ptsGauss(~id,:);
                       NSlave = NSlave(~id,:);
                       NSlaveMult = NSlaveMult(~id,:);
                    end
                 end
              end
           t = toc;
           M = M(obj.nodesSlave, obj.nodesMaster);
           D = D(obj.nodesSlave, obj.nodesSlave);
           if strcmp(mult_type,'dual')
              % make sure D is diagonal
              D = diag(sum(D,2));
           end
           
           if nargout == 3
              varargout{1} = t;
           elseif nargout == 4
              varargout{1} = t;
              varargout{2} = D\M;
           end
        end
   
        function getMatricesSize(obj)
            if isempty(obj.mshSlave)
                obj.nSMat = obj.nNodesSlave;
            else
                obj.nSMat = obj.mshSlave.nNodes;
            end

            if isempty(obj.mshMaster)
                obj.nMMat = obj.nNodesMaster;
            else
                obj.nMMat = obj.mshMaster.nNodes;
            end
        end


        function xiMaster = projectGP(obj,elM,elS,xi,n,xInt)
            xiMaster = zeros(length(xi),1);
            tol = 1e-8;
            itMax = 8;
            nodeM = obj.masterTopol(elM,:);
            nodeS = obj.slaveTopol(elS,:);
            coord = obj.masterCoord(:,1:2);
            % evaluate the basis functions on the integration point
            for ii = 1:length(xi)
                iter = 0;
                w = 0;
                Nslave = computeBasis1D(xi(ii),obj.degree);
                Nmaster = computeBasis1D(xiMaster(ii),obj.degree);
                % evaluate normal at integration point
                rhs1 = Nmaster*coord(nodeM,:);
                nS = Nslave*n(nodeS,:);
                rhs = rhs1' - w*nS' - xInt(ii,:)';
                % project Gauss Point
                while (norm(rhs,2) > tol) && (iter < itMax)
                    % compute Jacobian of the projection
                    iter = iter + 1;
                    J1 = computeJacProj(obj,nodeM,xi(ii));
                    J = [J1, -nS'];
                    ds = J\(-rhs);
                    xiMaster(ii) = xiMaster(ii) + ds(1);
                    w = w + ds(2);
                    Nmaster = computeBasis1D(xiMaster(ii),obj.degree);
                    % evaluate normal at integration point
                    rhs1 = Nmaster*coord(nodeM,:);
                    rhs = rhs1' - w*nS' - xInt(ii,:)';
                end

                if iter >= itMax
                    xiMaster(ii) = nan;
                    % mark gauss point that did not converge
                end
            end
        end
        %

        function J = computeJacProj(obj,nM,xi)
            % compute first colum of the Jacobian of the projection
            % derivative of basis functions on local unknown point xi
            % times coordinate of nodes
            c = obj.masterCoord(nM,1:2);
            switch obj.degree
                case 1
                    dN = [-0.5;0.5];
                    J = c'*dN;
                case 2
                    dN = [-0.5+xi;0.5+xi;-2*xi];
                    J = c'*dN;
            end
        end

        function Ndual = computeDualBasisF(obj,gpCoord)
            % compute dual basis function on slave GPs
            % get matrix of BF phi = Nslave*A
            % A = M\D
%             Nslave = computeBasis1D(gpPos,obj.degree);
%             M = Nslave'*(Nslave.*gpW);
%             M = M*(0.5*h);
%             D = diag(Nslave'*gpW)*(0.5*h);
%             A = M\D;
%             Nslave_dual = Nslave*A;   
             % use analytical definition from Popp et. al (only for planar
             % interfaces)
             switch obj.degree
               case 1
                 Ndual(:,1) = 0.5*(1-3*gpCoord);
                 Ndual(:,2) = 0.5*(1+3*gpCoord);
               case 2
                 N1 = @(x) 0.5*x.*(x-1);
                 N2 = @(x) 0.5*x.*(x+1);
                 N3 = @(x) (1-x).*(1+x);
                 Ndual(:,1) = N1(gpCoord) - 3/4*N3(gpCoord) + 0.5;
                 Ndual(:,2) = N2(gpCoord) -3/4*N3(gpCoord) + 0.5;
                 Ndual(:,3) = 5/2*N3(gpCoord) - 1;
             end
        end

        function n_n = computeNodalNormal(obj)
            % Return vector n of weighted normals
            topol = obj.slaveTopol;
            coord = obj.slaveCoord;
            % compute tangent of each element
            t = [coord(topol(:,2),1) - coord(topol(:,1),1), coord(topol(:,2),2) - coord(topol(:,1),2)];
            t = t./sqrt((t(:,1).^2 + t(:,2).^2));
            n = [t(:,2), -t(:,1)]; % chosen arbitrarily, later this will be checked

            % compute length of each element
            l = vecnorm([coord(topol(:,1),1) - coord(topol(:,2),1), coord(topol(:,1),2) - coord(topol(:,2),2)],2,2);
            n_n = zeros(length(coord),2);

            for i = 1:length(coord)
                el = find(any((ismember(topol,i)),2));
                if length(el) > 1
                    n_n(i,:) = (l(el(1))*n(el(1),:) + l(el(2))*n(el(2),:))/(l(el(1))+l(el(2)));
                elseif length(el) == 1
                    n_n(i,:) = n(el,:);
                elseif isempty(el)
                    n_n(i,:) = 0;
                end
            end
        end
        %

        function [wF,w1,pts] = getWeights(obj,nInt,type)
            switch obj.degree
                case 1
                    wF = zeros(nInt,obj.nElMaster*obj.nNelem);
                    w1 = zeros(nInt,obj.nElMaster);
                    pts = zeros(nInt,obj.nElMaster*2);
                    for i = 1:obj.nElMaster
                        [f, ptsInt] = computeMortarBasisF(obj,i, nInt, obj.masterTopol, obj.masterCoord);
                        fiMM = obj.computeRBFfiMM(ptsInt,type);
                        % solve local system to get weight of interpolant
                        wF(:,repNum(obj.nNelem,i)) = fiMM\f;
                        w1(:,i) = fiMM\ones(size(ptsInt,1),1);
                        pts(:,repNum(2,i)) = ptsInt;
                    end
                case 2
                    wF = zeros(nInt,obj.nElMaster*(obj.nNelem+1));
                    w1 = zeros(nInt,obj.nElMaster);
                    pts = zeros(nInt,obj.nElMaster*2);
                    for i = 1:obj.nElMaster
                        [f,ptsInt,fL] = computeMortarBasisF(obj,i, nInt, obj.masterTopol, obj.masterCoord);
                        fiMM = obj.computeRBFfiMM(ptsInt,type);
                        % solve local system to get weight of interpolant
                        wF(:,repNum(obj.nNelem+1,i)) = fiMM\[f,fL];
                        w1(:,i) = fiMM\ones(size(ptsInt,1),1);
                        pts(:,repNum(2,i)) = ptsInt;
                    end
            end
        end

        function [bf,pos,varargout] = computeMortarBasisF(obj,elemID,nInt,top,coord)
            intPts = [-1 1];
            intPts = linspace(intPts(1), intPts(2), nInt);
            % if strcmp(type,'wendland')
            %     intPts = sin(pi*intPts/2);
            % end
            bf = computeBasis1D(intPts,obj.degree);
            i1 = coord(top(elemID,1),:);
            i2 = coord(top(elemID,2),:);
            pos = ref2nod(intPts, i1, i2);
            if nargout == 3
                assert(obj.degree==2,'Incorrect number of outputs for basis functions of degree 1')
                bfL = computeBasis1D(intPts,1);
                varargout{1} = bfL(:,1);
            end
        end
        

        function [errNorm, fInt] = computeInterpError(obj,E,f)
           % f: function handle for analytical solution
            if E == 0
                errNorm = 0;
                fInt = 0;
                return
            end
            % compute quadratic error of interpolation on the slave side
            fM = f(obj.masterCoord(obj.nodesMaster,1)); % analytical function computed on master mesh
            fS = f(obj.slaveCoord(obj.nodesSlave,1)); % analytical function computed on master mesh
            lNod = computeLengthNodes(obj);
            % Quadratic error of interpolation for 1D mortar benchmarks
            fInt = E * fM;
            err2 = (fS - fInt).^2;
            errNorm = sqrt(sum(err2.*lNod(obj.nodesSlave)));
        end
        %
        function id = getWeightsID(obj,i)
            switch obj.degree
                case 1
                    id = [2*i-1 2*i];
                case 2
                    id = [3*i-2 3*i-1 3*i];
            end
        end

        function elemConnectivity = computeElementConnectivity(obj)
            % find connectivity of 1D interfaces
            % (trivially done looking at the x-coordinates only)
            elemConnectivity = zeros(obj.nElMaster,obj.nElSlave);
            for i = 1:obj.nElMaster
                tmp = sort([obj.masterCoord(obj.masterTopol(i,1),1),obj.masterCoord(obj.masterTopol(i,2),1)]);
                a = tmp(1); b = tmp(2);
                % loop trough master element to find connectivity
                for j = 1:obj.nElSlave
                    tmp = sort([obj.slaveCoord(obj.slaveTopol(j,1),1),obj.slaveCoord(obj.slaveTopol(j,2),1)]);
                    c = tmp(1); d = tmp(2);
                    if ~any([a>d,c>b])
                        % intersecting
                        elemConnectivity(i,j) = 1;
                    end
                end
            end
        end

        function lNod = computeLengthNodes(obj)
            lNod = zeros(obj.nNodesSlave,1);
            slave = obj.slaveCoord;
            switch obj.degree
                case 1
                    for el = obj.elemsSlave
                        % get length of the element
                        n1 = obj.slaveTopol(el,1);
                        n2 = obj.slaveTopol(el,2);
                        l = sqrt((slave(n1,1) - slave(n2,1))^2+(slave(n1,2) - slave(n2,2))^2);
                        lNod(n1) = lNod(n1) + l/2;
                        lNod(n2) = lNod(n2) + l/2;
                    end
                case 2
                    for el = 1:obj.nElSlave
                        % get length of the element
                        n1 = obj.slaveTopol(el,1);
                        n2 = obj.slaveTopol(el,2);
                        n3 = obj.slaveTopol(el,3);
                        l1 = sqrt((slave(n1,1) - slave(n2,1))^2+(slave(n1,2) - slave(n2,2))^2);
                        l2 = sqrt((slave(n2,1) - slave(n3,1))^2+(slave(n2,2) - slave(n3,2))^2);
                        lNod(n1) = lNod(n1) + l1/2;
                        lNod(n2) = lNod(n2) + l1/2 + l2/2;
                        lNod(n3) = lNod(n3) + l2/2;
                    end
            end
        end
        %



    end


    methods (Static)
        function r1 = computeRBFRadius(ptsInt)
            r1 = sqrt((max(ptsInt(:,1)) - min(ptsInt(:,1)))^2 + (max(ptsInt(:,2)) - min(ptsInt(:,2)))^2);
            r2 = (10/(length(ptsInt)-1))*r1;
            r = min(r1,r2);
        end

        function fiMM = computeRBFfiMM(ptsInt,type)
            r =  Mortar2D.computeRBFRadius(ptsInt);
            d = sqrt((ptsInt(:,1) - ptsInt(:,1)').^2 + (ptsInt(:,2) - ptsInt(:,2)').^2);
            fiMM = Mortar2D.rbfInterp(d,r,type);
        end

        function [fiNM] = computeRBFfiNM(ptsInt,ptsGauss,type)
            % id: id of points that has a distance < r with at least one
            % master points
            d = sqrt((ptsGauss(:,1) - ptsInt(:,1)').^2 + (ptsGauss(:,2) - ptsInt(:,2)').^2);
            r = sqrt((max(ptsInt(:,1)) - min(ptsInt(:,1)))^2 + (max(ptsInt(:,2)) - min(ptsInt(:,2)))^2);
            fiNM =  Mortar2D.rbfInterp(d,r,type);
        end

        function rbf = rbfInterp(d,r,type)
            % compute row of rbf interpolation matrix
            switch type
                case 'wendland'
                  v = 1-d./r;
                  v(v<0) = 0;
                    rbf = v.^4.*(1+d./r);
                case 'gauss'
                    rbf = exp(-d.^2/r^2);
                case 'imq'
                    rbf = (d.^2+r^2).^(-0.5);
            end
        end

    end
end

