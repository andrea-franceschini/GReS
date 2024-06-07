classdef Mortar3D < handle
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
        nodesMaster
        nodesSlave
        Dmat
        nSMat
        nMMat
        intMaster
        intSlave
    end

    methods
        function obj = Mortar3D(deg,varargin)
            % Constructor for the class
            % INPUT:
            % empty varargin (later manual assignment of
            % coordinates and topology of master and slave
            % non empty varargin: meshMaster, id master interface,
            % meshSlave, id slave interface
            obj.degree = deg;
            switch  nargin
                case 5
                    masterTop = varargin{1}.surfaces;
                    obj.masterCoord = varargin{1}.coordinates;
                    obj.masterTopol = masterTop(varargin{1}.surfaceTag == varargin{2}, :);
                    slaveTop = varargin{3}.surfaces;
                    obj.slaveCoord = varargin{3}.coordinates;
                    obj.slaveTopol = slaveTop(varargin{3}.surfaceTag == varargin{4}, :);
                    % get interfaces as mesh objects
                    obj.intMaster = varargin{1}.getSurfaceMesh(varargin{2});
                    obj.intSlave = varargin{3}.getSurfaceMesh(varargin{4});
                    getMatricesSize(obj,varargin{1},varargin{3});
                case 3 % cartgrid input
                    assert(varargin{1}.cartGrid && varargin{1}.nDim == 2,['Wrong number' ...
                        'of input for 2D Cartesian mesh']);
                    obj.masterCoord = varargin{1}.coordinates;
                    obj.masterTopol = varargin{1}.surfaces;
                    obj.slaveCoord = varargin{2}.coordinates;
                    obj.slaveTopol = varargin{2}.surfaces;
                    obj.intMaster = varargin{1};
                    obj.intSlave = varargin{2};
                    % elseif nargin == 6 && strcmp(varargin{1},'set') && obj.mshMaster.nDim == 2
                    %     obj.masterTopol = varargin{2};
                    %     obj.slaveTopol = varargin{3};
                    %     obj.masterCoord = varargin{4};
                    %     obj.slaveCoord = varargin{5};
                otherwise
                    error('Wrong number of inputs for class Mortar1D')
            end
            obj.nodesMaster = unique(obj.masterTopol);
            obj.nodesSlave = unique(obj.slaveTopol);
            obj.nElMaster = size(obj.masterTopol,1);
            obj.nElSlave = size(obj.slaveTopol,1);
            obj.nNodesMaster = size(obj.masterCoord,1);
            obj.nNodesSlave = size(obj.slaveCoord,1);
            if nargin == 3 % Set operator size for Cartesian meshes
                getMatricesSize(obj);
            end
            % Contact search algorithm to get connectivity between surface
            % meshes
            cs = ContactSearching(obj.intMaster,obj.intSlave,18);
            obj.elemConnectivity = cs.elemConnectivity;
            computeSlaveMatrix(obj);
        end

        function obj = computeSlaveMatrix(obj)
            % Gauss integration for slave mass matrix
            gaussMass = Gauss(12,3,2);
            elemSlave = Elements(obj.intSlave, gaussMass);
            l1 = 0;
            s1 = 16;
            for el = 1:size(obj.intSlave.surfaces,1)
                N1 = elemSlave.quad.getBasisFinGPoints();
                dJWeighed = elemSlave.quad.getDerBasisFAndDet(el,3);
                DLoc = N1'*diag(dJWeighed)*N1;
                dof = obj.slaveTopol(el,:);
                [jjLoc,iiLoc] = meshgrid(dof,dof);
                iiVec(l1+1:l1+s1) = iiLoc(:);
                jjVec(l1+1:l1+s1) = jjLoc(:);
                DVec(l1+1:l1+s1) = DLoc(:);
                l1 = l1 + s1;
            end
            % extract only mass matrix entries belonging to the interface
            obj.Dmat = sparse(iiVec, jjVec, DVec);
        end


        function [E,Mdetect,varargout] = computeMortarElementBased(obj,nGP)
            tic
            Mdetect = zeros(obj.nElMaster,obj.nElSlave);
            g = Gauss(12,nGP,2);
            %gpRef = g.coord;
            elemSlave = Elements(obj.intSlave,g);
            %gpWeight = g.weight;
            n = computeNodalNormal(obj,elemSlave);
            M = zeros(obj.nMMat,obj.nSMat);
            D = zeros(obj.nSMat,obj.nSMat);
            for i = 1:obj.nElSlave
                gpPos = g.coord; 
                NSlave = getBasisFinGPoints(elemSlave.quad); % Slave basis functions
                dJWeighed = elemSlave.quad.getDerBasisFAndDet(i,3); % Weighted Jacobian
                xSlave = getGPointsLocation(elemSlave.quad,i);
                idSlave = obj.slaveTopol(i,:);
                master_elems = find(obj.elemConnectivity(:,i)); 
                %id = true(length(gpRef),1);
                for m = master_elems'
                    [xiM] = obj.projectGP(m,i,gpPos,n,xSlave,elemSlave);
                    id = all([xiM >= -1, xiM <= 1],2);
                    Mdetect(m,i) = sum(id);
                    if any(id)
                        idMaster = obj.masterTopol(m,:);
                        NMaster = elemSlave.quad.computeBasisF(xiM(id,:));
                        if size(NMaster,2)~=4
                            NMaster = NMaster';
                            % this temporarily fix a bug generated by
                            % bsxfun when input is a scalar
                        end
                        Mloc = NMaster'*(NSlave(id,:).*dJWeighed(id)');
                        M(idMaster, idSlave) = M(idMaster, idSlave) + Mloc;
                        Dloc = NSlave(id,:)'*(NSlave(id,:).*dJWeighed(id)');
                        D(idSlave, idSlave) = D(idSlave, idSlave) + Dloc;
                        % sort out Points already projected
                        gpPos = gpPos(~id,:);
                        dJWeighed = dJWeighed(~id);
                        xSlave = xSlave(~id,:);
                        NSlave = NSlave(~id,:);
                    end
                end
            end
            t = toc;
            M = (M(obj.nodesMaster, obj.nodesSlave))';
            E = D\M;
            if nargout == 3
                varargout{1} = M;
            elseif nargout == 4
                varargout{1} = M;
                varargout{2} = t;
            end
        end

        function [E,Mdetect,varargout] = computeMortarRBF(obj,nGP,nInt,type)
            % this code work only for hexa 8 (quad 4 on the interface)
            Mdetect = zeros(obj.nElMaster,obj.nElSlave);
            tic
            % set Gauss class
            gM = Gauss(12,3,2); % gauss class for Master element interpolation
            g = Gauss(12,nGP,2); % gauss class for slave integration
            elemMaster = Elements(obj.intMaster, gM);
            elemSlave = Elements(obj.intSlave, g);
            M = zeros(obj.nMMat,obj.nSMat);
            D = zeros(obj.nSMat,obj.nSMat);
            % Perform interpolation on the master side (computing weights
            % and interpolation coordinates)
            [wFMat,w1Mat,ptsIntMat] = getWeights(obj,'master',nInt,elemMaster,type);
            [wFMatS,w1MatS,ptsIntMatS] = getWeights(obj,'slave',nInt,elemSlave,type);
            % Interpolation for support detection
            [wFSupp,w1Supp] = getSuppWeight(obj);
            % Loop trough slave elements
            for j = 1:obj.nElSlave
                % Compute Slave quantities
                %NSlave = getBasisFinGPoints(elemSlave.quad); % Slave basis functions
                dJWeighed = elemSlave.quad.getDerBasisFAndDet(j,3); % Weighted Jacobian
                % get Gauss Points position in the real space
                ptsGauss = getGPointsLocation(elemSlave.quad,j);
                idSlave = obj.slaveTopol(j,:);
                % compute slave basis function (still using radial basis
                % interpolation)
                ptsIntS = ptsIntMatS(:,[3*j-2 3*j-1 3*j]);
                fiNMS = obj.computeRBFfiNM(ptsIntS,ptsGauss,type);
                NSlave = (fiNMS*wFMatS(:,[4*j-3 4*j-2 4*j-1 4*j]))./(fiNMS*w1MatS(:,j));
                master_elems = find(obj.elemConnectivity(:,j));
                for jm = master_elems'
                    idMaster = obj.masterTopol(jm,:);
                    ptsInt = ptsIntMat(:,[3*jm-2 3*jm-1 3*jm]);              
                    fiNM = obj.computeRBFfiNM(ptsInt,ptsGauss,type);          
                    NMaster = (fiNM*wFMat(:,[4*jm-3 4*jm-2 4*jm-1 4*jm]))./(fiNM*w1Mat(:,jm));          
                    % automatically detect supports computing interpolant
                    % of planes
                    id = detectSupports(obj,jm,ptsGauss,wFSupp(:,[2*jm-1 2*jm]),w1Supp(:,jm));
                    Mdetect(jm,j) = sum(id);
                    if any(id)
                        NMaster = NMaster(id,:);
                        Mloc = NMaster'*(NSlave(id,:).*dJWeighed(id)');
                        Dloc = NSlave(id,:)'*(NSlave(id,:).*dJWeighed(id)');
                        M(idMaster, idSlave) = M(idMaster, idSlave) + Mloc;
                        D(idSlave, idSlave) = D(idSlave, idSlave) + Dloc;
                        % sort out Points already projected
                        dJWeighed = dJWeighed(~id);    
                        ptsGauss = ptsGauss(~id,:);
                        NSlave = NSlave(~id,:);
                    end
                end
            end
            t = toc;
            M = (M(obj.nodesMaster, obj.nodesSlave))';
            tic
            E = D\M;
            tSist = toc;
            %scale projection operator to recover partition of unity
            %E = E./sum(E,2);
            if nargout == 3
                varargout{1} = M;
            elseif nargout == 4
                varargout{1} = M;
                varargout{2} = t;
            elseif nargout == 5
                varargout{1} = M;
                varargout{2} = t;
                varargout{3} = tSist;
            end
        end

        function getMatricesSize(obj,varargin)
            if isempty(varargin)
                obj.nSMat = obj.nNodesSlave;
                obj.nMMat = obj.nNodesMaster;
            else
                obj.nMMat = varargin{1}.nNodes;
                obj.nSMat = varargin{2}.nNodes;
            end
        end

        function xiMaster = projectGP(obj,elM,elS,xi,n,xInt,elem)
            % xi: location of Gauss point in the reference space
            % xInt: location of GP in the physical space
            % elem: istance of elements class for Basis Function evaluation
            xiMaster = zeros(size(xi,1),2);
            tol = 1e-8;
            itMax = 8;
            nodeM = obj.masterTopol(elM,:);
            nodeS = obj.slaveTopol(elS,:);
            coord = obj.masterCoord;
            % loop trough each integration point to project
            for ii = 1:size(xi,1)
                iter = 0;
                w = 0;
                Nslave = elem.quad.computeBasisF(xi(ii,:));
                Nmaster = elem.quad.computeBasisF(xiMaster(ii,:));
                % evaluate normal at integration point
                rhs1 = coord(nodeM,:)'*Nmaster;
                nS = n(nodeS,:)'*Nslave;
                rhs = rhs1 - w*nS - xInt(ii,:)';

                % project Gauss Point
                while (norm(rhs,2) > tol) && (iter < itMax)
                    % compute Jacobian of the projection
                    iter = iter + 1;
                    % get local derivatives of basis functions
                    dN = elem.quad.computeDerBasisF(xiMaster(ii,:));
                    J1 = dN*coord(nodeM,:);
                    J = [J1', -nS];
                    ds = J\(-rhs);
                    xiMaster(ii,:) = xiMaster(ii,:) + ds(1:2)';
                    w = w + ds(3);
                    Nmaster = elem.quad.computeBasisF(xiMaster(ii,:));
                    % evaluate normal at integration point
                    rhs1 = coord(nodeM,:)'*Nmaster;
                    rhs = rhs1 - w*nS - xInt(ii,:)';
                end

                if iter >= itMax
                    xiMaster(ii) = nan;
                    % mark gauss point that did not converge
                end
            end
        end
        %

        function n_n = computeNodalNormal(obj,elem)
            % Return vector n of weighted normals
            n_n = zeros(obj.intSlave.nNodes,3);
            area = elem.quad.findAreaAndCentroid(1:obj.intSlave.nSurfaces); % area of each cell
            elem_normal = elem.quad.computeNormal(1:obj.intSlave.nSurfaces);
            topol = obj.slaveTopol;
            for i = 1:length(obj.slaveCoord)
                % get elements sharing node i
                elems = find(any(ismember(topol,i),2));
                n_n(i,:) = (area(elems)'*elem_normal(elems,:))/sum(area(elems));
            end
        end
        %

        function [errNorm, fInt] = computeInterpError(obj,E,f)
            if E == 0
                errNorm = 0;
                fInt = 0;
                return
            end
            % compute quadratic error of interpolation on the slave side
            fM = f(obj.masterCoord(:,1)); % analytical function computed on master mesh
            fS = f(obj.slaveCoord(:,1)); % analytical function computed on master mesh
            lNod = zeros(obj.nNodesSlave,1);
            slave = obj.slaveCoord;
            for el = 1:size(obj.nElSlave,1)
                % get length of the element
                n1 = obj.slaveTopol(el,1);
                n2 = obj.slaveTopol(el,2);
                l = sqrt((slave(n1,1) - slave(n2,1))^2+(slave(n1,2) - slave(n2,2))^2);
                lNod(n1) = lNod(n1) + l/2;
                lNod(n2) = lNod(n2) + l/2;
            end
            % Quadratic error of interpolation for 1D mortar benchmarks
            fInt = E * fM;
            err2 = (fS - fInt).^2;
            errNorm = sqrt(sum(err2.*lNod));
        end
    %
    %
        function fiMM = computeRBFfiMM(obj,ptsInt,type)
            fiMM = zeros(length(ptsInt), length(ptsInt));
            % compute RBF weight of interpolant
            % local radius for radial basis function interpolation
            r = obj.computeRBFradius(ptsInt);
            %r = sqrt((max(ptsInt(:,1)) - min(ptsInt(:,1)))^2 + (max(ptsInt(:,2)) - min(ptsInt(:,2)))^2 + (max(ptsInt(:,3)) - min(ptsInt(:,3)))^2);
            for ii = 1:length(ptsInt)
                d = sqrt((ptsInt(:,1) - ptsInt(ii,1)).^2 + ((ptsInt(:,2)) - ptsInt(ii,2)).^2 + ((ptsInt(:,3)) - ptsInt(ii,3)).^2);
                rbf = obj.rbfInterp(d,r,type);
                fiMM(ii,:) = rbf;
            end
        end

        function fiNM = computeRBFfiNM(obj,ptsInt, ptsGauss, type)
            fiNM = zeros(size(ptsGauss,1), size(ptsInt,1));
            r = sqrt((max(ptsInt(:,1)) - min(ptsInt(:,1)))^2 + (max(ptsInt(:,2)) - min(ptsInt(:,2)))^2 + (max(ptsInt(:,3)) - min(ptsInt(:,3)))^2);
            % compute RBF weight of interpolant
            % local radius for radial basis function interpolation
            for jj = 1:size(ptsGauss,1)
                d = sqrt((ptsInt(:,1) - ptsGauss(jj,1)).^2 + ((ptsInt(:,2)) - ptsGauss(jj,2)).^2 + ((ptsInt(:,3)) - ptsGauss(jj,3)).^2); %%%% FIXXXXXXXXXXXXXXXX
                rbf = obj.rbfInterp(d,r,type);
                fiNM(jj,:) = rbf;
            end
        end

        function [wF,w1] = getSuppWeight(obj)
            wF = zeros(4,obj.nElMaster*2); % weights for support interesection
            w1 = zeros(4,obj.nElMaster); % weight for rescaling of support interesection
            % loop trough master elements and interpolate Basis function
            for i = 1:obj.nElMaster
                ptsInt = obj.masterCoord(obj.masterTopol(i,:),:); % get real coordinates of master nodes
                fiMM = obj.computeRBFfiMM(ptsInt,'gauss');
                % solve local system to get weight of interpolant
                b = [1 1 0 0; 0 1 1 0];
                wF(:,[2*i-1 2*i]) = fiMM\b';
                w1(:,i) = fiMM\ones(size(ptsInt,1),1);
            end
        end

        function id = detectSupports(obj,iM,ptsG,wF,w1)
            ptsMaster = obj.masterCoord(obj.masterTopol(iM,:),:);
            fiNM = obj.computeRBFfiNM(ptsMaster,ptsG,'gauss');
            Nsupp = (fiNM*wF)./(fiNM*w1);
            id = all([Nsupp >= 0, Nsupp <= 1],2);
        end

        function [wF,w1,pts] = getWeights(obj,interface,nInt,elem,type)
            numPts = (nInt)^2;
            switch interface
                case 'master'
                    wF = zeros(numPts,obj.nElMaster*4);
                    w1 = zeros(numPts,obj.nElMaster);
                    pts = zeros(numPts,obj.nElMaster*3);
                    for i = 1:obj.nElMaster
                        [f, ptsInt] = computeBasisF2D(i, nInt, obj.masterTopol, obj.masterCoord, elem);
                        fiMM = obj.computeRBFfiMM(ptsInt,type);
                        % solve local system to get weight of interpolant
                        wF(:,[4*i-3 4*i-2 4*i-1 4*i]) = fiMM\f;
                        w1(:,i) = fiMM\ones(size(ptsInt,1),1);
                        pts(:,[3*i-2 3*i-1 3*i]) = ptsInt;
                    end
                case 'slave'
                    wF = zeros(numPts,obj.nElSlave*4);
                    w1 = zeros(numPts,obj.nElSlave);
                    pts = zeros(numPts,obj.nElSlave*3);
                    for i = 1:obj.nElSlave
                        [f, ptsInt] = computeBasisF2D(i, nInt, obj.slaveTopol, obj.slaveCoord, elem);
                        fiMM = obj.computeRBFfiMM(ptsInt,type);
                        % solve local system to get weight of interpolant
                        wF(:,[4*i-3 4*i-2 4*i-1 4*i]) = fiMM\f;
                        w1(:,i) = fiMM\ones(size(ptsInt,1),1);
                        pts(:,[3*i-2 3*i-1 3*i]) = ptsInt;
                    end
            end

            % loop trough master elements and interpolate Basis function

        end

    end

   


    methods(Static)
        function rbf = rbfInterp(d,r,type)
            % compute row of rbf interpolation matrix
            switch type
                case 'wendland'
                    rbf = pos(1-d./r).^4.*(1+d./r);
                case 'gauss'
                    rbf = exp(-d.^2/r^2);
                case 'imq'
                    rbf = (d.^2+r^2).^(-0.5);
            end
        end

        function r1 = computeRBFradius(ptsInt)
            % compute the radius of the RBF interpolation based on
            % coordinates of interpolation points
            r1 = sqrt((max(ptsInt(:,1)) - min(ptsInt(:,1)))^2 + (max(ptsInt(:,2)) - min(ptsInt(:,2)))^2 + (max(ptsInt(:,3)) - min(ptsInt(:,3)))^2);
            r2 = (5/sqrt(length(ptsInt)))*r1;
            r = min(r1,r2);
        end
    end
end

