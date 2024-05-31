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
        nodesMaster
        nodesSlave
        Dmat
        nSMat
        nMMat
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
            obj.nodesMaster = unique(obj.masterTopol);
            obj.nodesSlave = unique(obj.slaveTopol);
            obj.nElMaster = size(obj.masterTopol,1);
            obj.nElSlave = size(obj.slaveTopol,1);
            obj.nNodesMaster = size(obj.masterCoord,1);
            obj.nNodesSlave = size(obj.slaveCoord,1);
            obj.elemConnectivity = computeElementConnectivity(obj);
            getMatricesSize(obj);
            computeSlaveMatrix(obj);
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



        function [E, varargout] = computeMortarSegmentBased(obj,nGP)
            tic
            if length(unique(obj.slaveCoord(obj.nodesSlave,2))) ~= 1
                E = 0;
                return
            end
            g = Gauss(12,nGP,1);
            gpRef = g.coord;
            gpWeight = g.weight;
            M = zeros(obj.nMMat,obj.nSMat);
            for i = 1:obj.nElMaster
                % get element of the support
                a = obj.masterCoord(obj.masterTopol(i,1),:);
                b = obj.masterCoord(obj.masterTopol(i,2),:);
                % get shape function values and interpolation points coordinate
                slave_elems = find(obj.elemConnectivity(i,:));
                % loop trough connected slave elements
                for j = slave_elems
                    % get nodes
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
                    basisMaster = 0.5 + gMaster*[-0.5 0.5];
                    % compute basis function on slave elements points
                    gSlave = nod2ref(gPts, c, d);
                    basisSlave = 0.5 + gSlave.*[-0.5 0.5];
                    basisSlave = basisSlave.*gpWeight;
                    matVal = basisMaster'*basisSlave;
                    matVal = matVal.*(0.5*(i2(1)-i1(1)));
                    idMaster = [obj.masterTopol(i,1); obj.masterTopol(i,2)];
                    idSlave = [obj.slaveTopol(j,1); obj.slaveTopol(j,2)];
                    M(idMaster, idSlave) = M(idMaster, idSlave) + matVal;
                end
            end
            t = toc;
            M = (M(obj.nodesMaster, obj.nodesSlave))';
            E = obj.D\M;
            if nargout == 2
                varargout{1} = M;
            elseif nargout == 3
                varargout{1} = M;
                varargout{2} = t;
            end
        end

        function [E, varargout] = computeMortarElementBased(obj,nGP)
            tic
            g = Gauss(12,nGP,1);
            gpRef = g.coord;
            gpWeight = g.weight;
            n = computeNodalNormal(obj);
            M = zeros(obj.nMMat,obj.nSMat);
            D = zeros(obj.nSMat,obj.nSMat);
            nN = 2; % numb. nodes per element
            if obj.degree > 1
                nN = 3;
            end
            for i = 1:obj.nElSlave
                % get nodes
                gpPos = gpRef;
                gpW = gpWeight;
                s1 = obj.slaveCoord(obj.slaveTopol(i,1),:);
                s2 = obj.slaveCoord(obj.slaveTopol(i,nN),:);
                idSlave = obj.slaveTopol(i,1:nN);
                xSlave = ref2nod(gpPos, s1, s2);
                master_elems = find(obj.elemConnectivity(:,i));
                %id = true(length(gpRef),1);
                for m = master_elems'
                    [xiM] = obj.projectGP(m,i,gpPos,n,xSlave);
                    id = all([xiM >= -1, xiM <= 1],2);
                    if any(id)
                        NMaster = computeBasis(obj,xiM(id));
                        NSlave = computeBasis(obj,gpPos(id));
                        Mloc = NMaster'*(NSlave.*gpW(id));
                        Mloc = Mloc*(0.5*sqrt((s1(1)-s2(1))^2 + (s1(2)-s2(2))^2));
                        Dloc = NSlave'*(NSlave.*gpW(id));
                        Dloc = Dloc*(0.5*sqrt((s1(1)-s2(1))^2 + (s1(2)-s2(2))^2));
                        idMaster = obj.masterTopol(m,:);
                        M(idMaster, idSlave) = M(idMaster, idSlave) + Mloc;
                        D(idSlave, idSlave) = D(idSlave, idSlave) + Dloc;
                        % sort out Points already projected
                        gpPos = gpPos(~id);
                        gpW = gpW(~id);
                        xSlave = xSlave(~id,:);
                    end
                end
            end
            t = toc;
            M = (M(obj.nodesMaster, obj.nodesSlave))';
            E = D\M;
            if nargout == 2
                varargout{1} = M;
            elseif nargout == 3
                varargout{1} = M;
                varargout{2} = t;
            end
        end

        function [E,varargout] = computeMortarRBF(obj,nGP,nInt,type)
            % type: family of Radial Basis function to use
            tic
            g = Gauss(12,nGP,1);
            gpRef = g.coord;
            gpWeight = g.weight;
            M = zeros(obj.nMMat,obj.nSMat);
            D = zeros(obj.nSMat,obj.nSMat);
            n = nInt;
            nN = 2; % numb. nodes per element
            if obj.degree > 1
                n = 2*nInt-1;
                nN = 3;
            end
            wFMat = zeros(n,obj.nElMaster*nN);
            w1Mat = zeros(n,obj.nElMaster);
            ptsIntMat = zeros(n,obj.nElMaster*2);
            if obj.degree > 1
                wFaux = zeros(4,obj.nElMaster);
                w1aux = zeros(4,obj.nElMaster);
                ptsaux = zeros(4,obj.nElMaster*2);
            end
            for i = 1:obj.nElMaster
                % get shape function values and interpolation points
                % coordinates
                [vals, ptsInt] = computeBasisF1D(i, nInt, obj.masterTopol, obj.masterCoord, obj.degree);
                % auxiliary interpolation for contact detection
                if obj.degree > 1
                    [v,p] = computeBasisF1D(i, nInt, obj.masterTopol, obj.masterCoord, 0);
                    fiMMaux = computeRBFfiMM(obj,p,type);
                    wFaux(:,i) = fiMMaux\v;
                    w1aux(:,i) = fiMMaux\ones(size(p,1),1);
                    ptsaux(:,[2*i-1 2*i]) = p;
                end
                fiMM = computeRBFfiMM(obj,ptsInt,type);
                % solve local system to get weight of interpolant. Store
                % interpolation results into matrices
                wFMat(:,getWeightsID(obj,i)) = fiMM\vals;
                w1Mat(:,i) = fiMM\ones(size(ptsInt,1),1);
                ptsIntMat(:,[2*i-1 2*i]) = ptsInt;
            end
            % Loop trough slave elements
            for j = 1:obj.nElSlave
                % Gauss parameters
                gpPos = gpRef;
                gpW = gpWeight;
                % nodes of slave element
                s1 = obj.slaveCoord(obj.slaveTopol(j,1),:);
                s2 = obj.slaveCoord(obj.slaveTopol(j,nN),:);
                idSlave = obj.slaveTopol(j,1:nN);
                % Real position of Gauss points
                ptsGauss = ref2nod(gpPos, s1, s2);
                % Get connected master elements
                master_elems = find(obj.elemConnectivity(:,j));
                % Loop troug master elements
                for jm = master_elems'
                    ptsInt = ptsIntMat(:,[2*jm-1 2*jm]);
                    idMaster = obj.masterTopol(jm,1:nN);
                    fiNM = computeRBFfiNM(obj,ptsInt,ptsGauss,type);
                    %NMaster = (fiNM*wFMat(:,getWeightsID(obj,jm)))./(fiNM*w1Mat(:,jm));
                    NMaster = fiNM*wFMat(:,getWeightsID(obj,jm));
                    if obj.degree > 1
                        % do the contact checking for quadratic elements
                        p = ptsaux(:,[2*jm-1 2*jm]);
                        fiNMaux = computeRBFfiNM(obj,p,ptsGauss,type);
                        Naux = (fiNMaux*wFaux(:,jm))./(fiNMaux*w1aux(:,jm));
                        id = all([Naux >= 0, Naux <= 1],2); % contact detection
                    else
                        id = all([NMaster >= 0, NMaster <= 1],2); % contact detection
                    end
                    if any(id)
                        NMaster = NMaster(id,:);
                        NSlave = computeBasis(obj,gpPos(id));
                        Mloc = NMaster'*(NSlave.*gpW(id));
                        Mloc = Mloc*(0.5*sqrt((s1(1)-s2(1))^2 + (s1(2)-s2(2))^2));
                        Dloc = NSlave'*(NSlave.*gpW(id));
                        Dloc = Dloc*(0.5*sqrt((s1(1)-s2(1))^2 + (s1(2)-s2(2))^2));
                        M(idMaster, idSlave) = M(idMaster, idSlave) + Mloc;
                        D(idSlave, idSlave) = D(idSlave, idSlave) + Dloc;
                        % sort out Points already projected
                        gpPos = gpPos(~id);    % errNormEB(sizeCount) = computeInterpError(E_EB,fMaster,fSlave,lNod);~id);
                        gpW = gpW(~id);
                        ptsGauss = ptsGauss(~id,:);
                    end
                end
            end
            t = toc;
            M = (M(obj.nodesMaster, obj.nodesSlave))';
            E = D\M;
            if nargout == 2
                varargout{1} = M;
            elseif nargout == 3
                varargout{1} = M;
                varargout{2} = t;
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

        function N = computeBasis(obj,pos)
            if obj.degree > 1
                N = computeQuadraticSF(pos);
            else
                N = compute1DBasisF(pos);
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
                Nslave = computeBasis(obj,xi(ii));
                Nmaster = computeBasis(obj,xiMaster(ii));
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
                    Nmaster = computeBasis(obj,xiMaster(ii));
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
            c = obj.masterCoord(nM,:);
            switch obj.degree
                case 1
                    dN = [-0.5;0.5];
                    J = c'*dN;
                case 2
                    dN = [-0.5+xi;-2*xi;0.5+xi];
                    J = c'*dN;
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

        function [errNorm, fInt] = computeInterpError(obj,E,f)
            if E == 0
                errNorm = 0;
                fInt = 0;
                return
            end
            % compute quadratic error of interpolation on the slave side
            fM = f(obj.masterCoord(:,1)); % analytical function computed on master mesh
            fS = f(obj.slaveCoord(:,1)); % analytical function computed on master mesh
            lNod = computeLengthNodes(obj);
            % Quadratic error of interpolation for 1D mortar benchmarks
            fInt = E * fM;
            err2 = (fS - fInt).^2;
            errNorm = sqrt(sum(err2.*lNod));
        end
        %

        function rbf = rbfInterp(obj,d,r,type)
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

        function r1 = computeRBFRadius(obj,ptsInt)
            r1 = sqrt((max(ptsInt(:,1)) - min(ptsInt(:,1)))^2 + (max(ptsInt(:,2)) - min(ptsInt(:,2)))^2);
            r2 = (5/(length(ptsInt)-1))*r1;
            r = min(r1,r2);
        end

        function fiMM = computeRBFfiMM(obj,ptsInt,type)
            % compute interpolation matrix
            fiMM = zeros(length(ptsInt), length(ptsInt));
            % compute RBF weight of interpolant
            % radius of RBF interpolation
            r = computeRBFRadius(obj,ptsInt);
            for ii = 1:length(ptsInt)
                d = sqrt((ptsInt(:,1) - ptsInt(ii,1)).^2 + ((ptsInt(:,2)) - ptsInt(ii,2)).^2);
                rbf = obj.rbfInterp(d,r,type);
                fiMM(ii,:) = rbf;
            end
        end

        function fiNM = computeRBFfiNM(obj,ptsInt, ptsGauss, type)
            % compute interpolation matrix
            fiNM = zeros(size(ptsGauss,1), size(ptsInt,1));
            % compute interpolant on local Gauss points
            r = computeRBFRadius(obj,ptsInt);
            for jj = 1:size(ptsGauss,1)
                d = sqrt((ptsInt(:,1) - ptsGauss(jj,1)).^2 + ((ptsInt(:,2)) - ptsGauss(jj,2)).^2);
                rbf = obj.rbfInterp(d,r,type);
                fiNM(jj,:) = rbf;
            end
        end

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
                    for el = 1:obj.nElSlave
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
end

