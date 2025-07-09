classdef Mortar3D < handle
    % Class implementing algorithms for mortar method with 1D interfaces 
    properties
        mshMaster
        mshSlave
        masterTopol
        slaveTopol
        masterCoord
        slaveCoord            
        nElMaster             % number of master elements
        nElSlave              % number of slave elements
        nNodesMaster          % number of master nodes
        nNodesSlave           % number of slave nodes
        elemConnectivity      % connectivity matrix
        degree                % degree of interpolation of mortar 
        nodesMaster           % master side node subset
        nodesSlave            % slave side node subset
        elemNormal
        %Dmat
        nSMat
        nMMat
        intMaster             % 2D mesh object refering to the master side
        intSlave              % 2D mesh object refering to the slave side
        nNmaster              % number of nodes per elements in master mesh
        nNslave               % number of nodes per elements in master mesh
        masterCellType
        slaveCellType
        D
        M
        wF                    % matrix of RBF weights
        w1                    % matrix of RBF weight for scaling
        ptsRBF                % RBF interpolation points
        edge2nodes
        edge2cells
        nEdges
        f2cMaster
        f2cSlave
        areaMap
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
                obj.mshMaster = varargin{1};
                obj.mshSlave = varargin{3};
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
          if deg == 1
             switch obj.intMaster.surfaceVTKType(1)
                case 5 % Triangle mesh Master
                   obj.nNmaster = 3;
                   obj.masterCellType = 10;
                case 9 % Quad mesh Master
                   obj.nNmaster = 4;
                   obj.masterCellType = 12;
             end
          else
             assert(obj.intMaster.surfaceVTKType(1)==9,['Quadratic triangle not' ...
                'yet implemented in Mortar class'])
             obj.nNmaster = 8;
          end
          switch obj.intSlave.surfaceVTKType(1)
             case 5 % Triangle mesh Master
                obj.slaveCellType = 10;
                obj.nNslave = 3;
             case 9 % Quad mesh Master
                obj.slaveCellType = 12;
                obj.nNslave = 4;
          end
          obj.nElMaster = size(obj.masterTopol,1);
          obj.nElSlave = size(obj.slaveTopol,1);
          obj.nNodesMaster = size(obj.masterCoord,1);
          obj.nNodesSlave = size(obj.slaveCoord,1);
          if nargin == 3 % Set operator size for Cartesian meshes
             getMatricesSize(obj);
          end
          % Contact search algorithm to get connectivity between surface
          % meshes
          if nargin == 3 && length(unique(obj.masterCoord(:,3)))<2
             if obj.degree < 2
                mshM = varargin{1};
                mshS = varargin{2};
             else
                mshM = getQuad4mesh(varargin{1});
                mshS = getQuad4mesh(varargin{2});
             end
             obj.elemConnectivity = cartConn(obj,mshM,mshS);
          else
             if length(unique(obj.masterCoord(:,3)))<2
                kdop = 8;
             else
                kdop = 18;
             end
             switch obj.degree
                case 1
                   cs = ContactSearching(obj.intMaster,obj.intSlave,kdop);
                case 2 % perform contact searching based on a quad4 mesh
                   intM = obj.intMaster.getQuad4mesh();
                   intS = obj.intSlave.getQuad4mesh();
                   cs = ContactSearching(intM,intS,kdop);
             end
             obj.elemConnectivity = cs.elemConnectivity;
          end
          % connected master and slave surfaces (find non empty rows and
          % columns of connectivity matrix)
          idM = sum(obj.elemConnectivity,2) > 0;
          idS = sum(obj.elemConnectivity,1) > 0;
          % get list of nodes that actually belong to elements in contact
          obj.nodesMaster = unique(obj.masterTopol(idM,:));
          obj.nodesSlave = unique(obj.slaveTopol(idS,:));
          setupEdgeTopology(obj);
          [obj.f2cMaster,obj.f2cSlave] = obj.buildFace2CellMap();
          %
          %computeSlaveMatrix(obj);
       end

       function obj = computeSlaveMatrix(obj)
          % Gauss integration for slave mass matrix
          gaussMass = Gauss(12,3,2);
          elemSlave = Elements(obj.intSlave, gaussMass);
          l1 = 0;
          s1 = size(obj.slaveTopol,2)^2;
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


       function [D,M,varargout] = computeMortarElementBased(obj,nGP)
          tic
          c_ns = 0;
          Mdetect = zeros(obj.nElMaster,obj.nElSlave);
          g = Gauss(obj.slaveCellType,nGP,2);
          %gpRef = g.coord;
          elemSlave = getElem(obj,g,'slave');
          elemMaster = getElem(obj,g,'master');
          %gpWeight = g.weight;
          n = computeNodalNormal(obj,elemSlave);
          [imVec,jmVec,Mvec] = deal(zeros(nnz(obj.elemConnectivity)*obj.nNmaster^2,1));
          [isVec,jsVec,Dvec] = deal(zeros(nnz(obj.elemConnectivity)*obj.nNslave^2,1));
          cs = 0;
          cm = 0;
          for i = 1:obj.nElSlave
             gpPos = g.coord;
             NSlave = getBasisFinGPoints(elemSlave); % Slave basis functions
             dJWeighed = elemSlave.getDerBasisFAndDet(i,3); % Weighted Jacobian
             xSlave = getGPointsLocation(elemSlave,i);
             idSlave = obj.slaveTopol(i,:);
             master_elems = find(obj.elemConnectivity(:,i));
             %id = true(length(gpRef),1);
             for m = master_elems'
                [xiM] = obj.projectGP(m,i,gpPos,n,xSlave,elemSlave,elemMaster);
                id = detectSupportEB(obj,xiM);
                Mdetect(m,i) = sum(id);
                if any(id)
                   idMaster = obj.masterTopol(m,:);
                   NMaster = elemMaster.computeBasisF(xiM(id,:));
                   Mloc = NSlave(id,:)'*(NMaster.*dJWeighed(id)');
                   Dloc = NSlave(id,:)'*(NSlave(id,:).*dJWeighed(id)');
                   nm = numel(Mloc);
                   ns = numel(Dloc);
                   % keeping M and D sparse to improve performance
                   [jjM,iiM] = meshgrid(idMaster,idSlave);
                   [jjS,iiS] = meshgrid(idSlave,idSlave);
                   imVec(cm+1:cm+nm) = iiM(:); 
                   jmVec(cm+1:cm+nm) = jjM(:);
                   isVec(cs+1:cs+ns) = iiS(:); 
                   jsVec(cs+1:cs+ns) = jjS(:);
                   Mvec(cm+1:cm+nm) = Mloc(:);
                   Dvec(cs+1:cs+ns) = Dloc(:);
                   % sort out Points already projected
                   gpPos = gpPos(~id,:);
                   dJWeighed = dJWeighed(~id);
                   xSlave = xSlave(~id,:);
                   NSlave = NSlave(~id,:);
                   cs = cs+ns;
                   cm = cm+nm;
                end
             end
             if ~all(id)
                fprintf('GP not sorted for slave elem %i \n',i);
                c_ns = c_ns + 1;
             end
          end
          t = toc;
            imVec = imVec(1:cm); jmVec = jmVec(1:cm);
            isVec = isVec(1:cs); jsVec = jsVec(1:cs);
            Mvec = Mvec(1:cm); Dvec = Dvec(1:cs);
            M = sparse(imVec,jmVec,Mvec,obj.nNodesSlave,obj.nNodesMaster);
            M = M(obj.nodesSlave, obj.nodesMaster);
            D = sparse(isVec,jsVec,Dvec,obj.nNodesSlave,obj.nNodesSlave);
            D = D(obj.nodesSlave,obj.nodesSlave);
            %E = D\M;
            if nargout == 3
                varargout{1} = t;
            elseif nargout == 4
                varargout{1} = t;
                varargout{2} = Mdetect;
            elseif nargout == 5
               varargout{1} = t;
                varargout{2} = Mdetect;
                varargout{3} = D\M;
            end
            fprintf('Points not sorted out: %i \n',c_ns)
       end

       function [D,M,varargout] = computeMortarRBF(obj,nGP,nInt,type,mult_type)
            c_ns = 0;
            Mdetect = zeros(obj.nElMaster,obj.nElSlave);
            % set Gauss class
            gM = Gauss(obj.masterCellType,3,2); % gauss class for Master element interpolation
            gS = Gauss(obj.slaveCellType,nGP,2); % gauss class for slave integration
            elemMaster = getElem(obj,gM,'master');
            elemSlave = getElem(obj,gS,'slave');
            [imVec,jmVec,Mvec] = deal(zeros(nnz(obj.elemConnectivity)*obj.nNmaster^2,1));
            [isVec,jsVec,Dvec] = deal(zeros(nnz(obj.elemConnectivity)*obj.nNmaster^2,1));
            % Perform interpolation on the master side (computing weights
            % and interpolation coordinates)
            getWeights(obj,'master',nInt,elemMaster,type);
            %[wFMatS,w1MatS,ptsIntMatS] = getWeights(obj,'slave',nInt,elemSlave,type);
            % Interpolation for support detection
            %[wFSupp,w1Supp] = getSuppWeight(obj);
            %tInterp = toc;
            % Loop trough slave elements
            cs = 0;
            cm = 0;
            for j = 1:obj.nElSlave
               %Compute Slave quantities
               dJWeighed = elemSlave.getDerBasisFAndDet(j,3); % Weighted Jacobian
               %get Gauss Points position in the real space
               ptsGauss = getGPointsLocation(elemSlave,j);
               idSlave = obj.slaveTopol(j,:);
               NSlave = getBasisFinGPoints(elemSlave); % Slave basis functions
               switch mult_type
                  case 'standard'
                     NSlaveMult = NSlave; % Slave basis functions
                  case 'dual'
                     NSlaveMult = obj.computeDualBasisF(NSlave,dJWeighed);
               end
               % compute slave basis function (still using radial basis
               % interpolation)
               %ptsIntS = ptsIntMatS(:,[3*j-2 3*j-1 3*j]);
               %fiNMS = obj.computeRBFfiNM(ptsIntS,ptsGauss,type);
               %NSlave = (fiNMS*wFMatS(:,[4*j-3 4*j-2 4*j-1 4*j]))./(fiNMS*w1MatS(:,j));
               master_elems = find(obj.elemConnectivity(:,j));
               for jm = master_elems'
                  idMaster = obj.masterTopol(jm,:);
                  [NMaster,id] = obj.getMasterBasis(jm,ptsGauss); % compute interpolated master basis function \Pi(Nm)
                  % ptsInt = ptsIntMat(:,repNum(3,jm));
                  % [fiNM,id1] = obj.computeRBFfiNM(ptsInt,ptsGauss,type);
                  % switch obj.degree
                  %     case 1
                  %         NMaster = (fiNM*wFMat(:,repNum(obj.nNmaster,jm)))./(fiNM*w1Mat(:,jm));
                  %         Nsupp = NMaster(:,[1 2 3]);
                  %     case 2
                  %         Ntmp = (fiNM*wFMat(:,repNum(obj.nNmaster+2,jm)))./(fiNM*w1Mat(:,jm));
                  %         NMaster = Ntmp(:,1:obj.nNmaster);
                  %         Nsupp = Ntmp(:,[end-1 end]);
                  % end
                  % % automatically detect supports computing interpolant
                  % id = all([Nsupp >= 0-tol id1],2);
                  Mdetect(jm,j) = sum(id);
                  if any(id)
                     if ~strcmp(mult_type,'P0')
                        NMaster = NMaster(id,:);
                        Mloc = NSlaveMult(id,:)'*(NMaster.*dJWeighed(id)');
                        Dloc = NSlaveMult(id,:)'*(NSlave(id,:).*dJWeighed(id)');
                        nm = numel(Mloc);
                        ns = numel(Dloc);
                        % keeping M and D sparse to improve performance
                        [jjM,iiM] = meshgrid(idMaster,idSlave);
                        [jjS,iiS] = meshgrid(idSlave,idSlave);
                        imVec(cm+1:cm+nm) = iiM(:); jmVec(cm+1:cm+nm) = jjM(:);
                        isVec(cs+1:cs+ns) = iiS(:); jsVec(cs+1:cs+ns) = jjS(:);
                        Mvec(cm+1:cm+nm) = Mloc(:);
                        Dvec(cs+1:cs+ns) = Dloc(:);
                        % sort out Points already projected
                        dJWeighed = dJWeighed(~id);
                        ptsGauss = ptsGauss(~id,:);
                        NSlave = NSlave(~id,:);
                        NSlaveMult = NSlaveMult(~id,:);
                        cs = cs+ns;
                        cm = cm+nm;
                     else % P0 MULTIPLIER
                        NMaster = NMaster(id,:);
                        Mloc = sum(NMaster.*dJWeighed(id)',1);
                        Dloc = sum(NSlave(id,:).*dJWeighed(id)',1);
                        nm = numel(Mloc);
                        ns = numel(Dloc);
                        % keeping M and D sparse to improve performance
                        [jjM,iiM] = meshgrid(idMaster,j);
                        [jjS,iiS] = meshgrid(idSlave,j);
                        imVec(cm+1:cm+nm) = iiM(:); jmVec(cm+1:cm+nm) = jjM(:);
                        isVec(cs+1:cs+ns) = iiS(:); jsVec(cs+1:cs+ns) = jjS(:);
                        Mvec(cm+1:cm+nm) = Mloc(:);
                        Dvec(cs+1:cs+ns) = Dloc(:);
                        % sort out Points already projected
                        dJWeighed = dJWeighed(~id);
                        ptsGauss = ptsGauss(~id,:);
                        NSlave = NSlave(~id,:);
                        cs = cs+ns;
                        cm = cm+nm;
                     end
                  end
               end
               if ~all(id)
                  fprintf('GP not sorted for slave elem %i \n',j);
                  c_ns = c_ns + 1;
               end
            end
            imVec = imVec(1:cm); jmVec = jmVec(1:cm);
            isVec = isVec(1:cs); jsVec = jsVec(1:cs);
            Mvec = Mvec(1:cm); Dvec = Dvec(1:cs);
            if strcmp(mult_type,'P0')
               M = sparse(imVec,jmVec,Mvec,obj.nElSlave,obj.nNodesMaster);
               M = M(:,obj.nodesMaster);
               D = sparse(isVec,jsVec,Dvec,obj.nElSlave,obj.nNodesSlave);
               D = D(:,obj.nodesSlave);
            else
               M = sparse(imVec,jmVec,Mvec,obj.nNodesSlave,obj.nNodesMaster);
               M = M(obj.nodesSlave, obj.nodesMaster);
               D = sparse(isVec,jsVec,Dvec,obj.nNodesSlave,obj.nNodesSlave);
               D = D(obj.nodesSlave,obj.nodesSlave);
            end
            if strcmp(mult_type,'dual')
               D = diag(sum(D,2)); % make sure D is diagonal by lumping
            end
            obj.D = D;
            obj.M = M;
            %scale projection operator to recover partition of unity
            %E = E./sum(E,2);
            if nargout == 3
                varargout{1} = [tInterp,tInteg];
            elseif nargout == 4
                varargout{1} = [tInterp,tInteg];
                varargout{2} = Mdetect;
            elseif nargout == 5
                tic
                varargout{2} = Mdetect;
                varargout{3} = D\M;
                tSist = toc;
                varargout{1} = [tInterp,tInteg,tSist];
            end
            fprintf('Points not sorted out: %i \n',c_ns)
        end

%         function [D,M,varargout] = computeMortarRBF(obj,nGP,nInt,type)
%             % this code work only for hexa 8 (quad 4 on the interface)
%             Mdetect = zeros(obj.nElMaster,obj.nElSlave);
%             tic
%             % set Gauss class
%             gM = Gauss(obj.masterCellType,3,2); % gauss class for Master element interpolation
%             gS = Gauss(obj.slaveCellType,nGP,2); % gauss class for slave integration
%             elemMaster = getElem(obj,gM,'master');
%             elemSlave = getElem(obj,gS,'slave');
%             [imVec,jmVec,Mvec] = deal(zeros(nnz(obj.elemConnectivity)*obj.nNmaster^2,1));
%             [isVec,jsVec,Dvec] = deal(zeros(nnz(obj.elemConnectivity)*obj.nNmaster^2,1));
%             % Perform interpolation on the master side (computing weights
%             % and interpolation coordinates)
%             tic
%             [wFMat,w1Mat,ptsIntMat] = getWeights(obj,'master',nInt,elemMaster,type);
%             %[wFMatS,w1MatS,ptsIntMatS] = getWeights(obj,'slave',nInt,elemSlave,type);
%             % Interpolation for support detection
%             %[wFSupp,w1Supp] = getSuppWeight(obj);
%             tInterp = toc;
%             tic
%             % Loop trough slave elements
%             cs = 0;
%             cm = 0;
%             for j = 1:obj.nElSlave
%                 % Compute Slave quantities
%                 dJWeighed = elemSlave.getDerBasisFAndDet(j,3); % Weighted Jacobian
%                 % get Gauss Points position in the real space
%                 ptsGauss = getGPointsLocation(elemSlave,j);
%                 idSlave = obj.slaveTopol(j,:);
%                 NSlave = getBasisFinGPoints(elemSlave); % Slave basis functions
% %                 switch mult_type
% %                     case 'standard'
% %                         NSlaveMult = NSlave; % Slave basis functions
% %                     case 'dual'
% %                         NSlaveMult = computeDualBasisF(obj,elemSlave,ptsGauss,dJWeighed);
% %                 end
%                 % compute slave basis function (still using radial basis
%                 % interpolation)
%                 %ptsIntS = ptsIntMatS(:,[3*j-2 3*j-1 3*j]);
%                 %fiNMS = obj.computeRBFfiNM(ptsIntS,ptsGauss,type);
%                 %NSlave = (fiNMS*wFMatS(:,[4*j-3 4*j-2 4*j-1 4*j]))./(fiNMS*w1MatS(:,j));
%                 master_elems = find(obj.elemConnectivity(:,j));
%                 for jm = master_elems'
%                     idMaster = obj.masterTopol(jm,:);
%                     ptsInt = ptsIntMat(:,repNum(3,jm));              
%                     [fiNM,id1] = obj.computeRBFfiNM(ptsInt,ptsGauss,type);
%                     switch obj.degree
%                         case 1
%                             NMaster = (fiNM*wFMat(:,repNum(obj.nNmaster,jm)))./(fiNM*w1Mat(:,jm));
%                             Nsupp = NMaster(:,[1 2 3]);
%                         case 2
%                             Ntmp = (fiNM*wFMat(:,repNum(obj.nNmaster+2,jm)))./(fiNM*w1Mat(:,jm));
%                             NMaster = Ntmp(:,1:obj.nNmaster);
%                             Nsupp = Ntmp(:,[end-1 end]);
%                     end
%                     % automatically detect supports computing interpolant
%                     id = all([Nsupp >= 0, Nsupp <= 1 id1],2);
%                     Mdetect(jm,j) = sum(id);
%                     if any(id)
%                         NMaster = NMaster(id,:);
%                         Mloc = NSlave(id,:)'*(NMaster.*dJWeighed(id)');
%                         Dloc = NSlave(id,:)'*(NSlave(id,:).*dJWeighed(id)');
%                         nm = numel(Mloc);
%                         ns = numel(Dloc);
%                         % keeping M and D sparse to improve performance
%                         [jjM,iiM] = meshgrid(idMaster,idSlave);
%                         [jjS,iiS] = meshgrid(idSlave,idSlave);
%                         imVec(cm+1:cm+nm) = iiM(:); jmVec(cm+1:cm+nm) = jjM(:);
%                         isVec(cs+1:cs+ns) = iiS(:); jsVec(cs+1:cs+ns) = jjS(:);
%                         Mvec(cm+1:cm+nm) = Mloc(:);
%                         Dvec(cs+1:cs+ns) = Dloc(:);
%                         % sort out Points already projected
%                         dJWeighed = dJWeighed(~id);    
%                         ptsGauss = ptsGauss(~id,:);
%                         NSlave = NSlave(~id,:);
%                         cs = cs+ns;
%                         cm = cm+nm;
%                     end
%                 end
%             end
%             tInteg = toc;
%             imVec = imVec(1:cm); jmVec = jmVec(1:cm);
%             isVec = isVec(1:cs); jsVec = jsVec(1:cs);
%             Mvec = Mvec(1:cm); Dvec = Dvec(1:cs);
%             M = sparse(imVec,jmVec,Mvec,obj.nNodesSlave,obj.nNodesMaster);
%             M = M(obj.nodesSlave, obj.nodesMaster);
%             D = sparse(isVec,jsVec,Dvec,obj.nNodesSlave,obj.nNodesSlave);
%             D = D(obj.nodesSlave,obj.nodesSlave);
%             %scale projection operator to recover partition of unity
%             %E = E./sum(E,2);
%             if nargout == 3
%                 varargout{1} = [tInterp,tInteg];
%             elseif nargout == 4
%                 varargout{1} = [tInterp,tInteg];
%                 varargout{2} = Mdetect;
%             elseif nargout == 5
%                 tic
%                 varargout{2} = Mdetect;
%                 varargout{3} = D\M;
%                 tSist = toc;
%                 varargout{1} = [tInterp,tInteg,tSist];
%             end
%         end

        function getMatricesSize(obj,varargin)
            if isempty(varargin)
                obj.nSMat = obj.nNodesSlave;
                obj.nMMat = obj.nNodesMaster;
            else
                obj.nMMat = varargin{1}.nNodes;
                obj.nSMat = varargin{2}.nNodes;
            end
        end

        function xiMaster = projectGP(obj,elM,elS,xi,n,xInt,elemSlave,elemMaster)
            % xi: location of Gauss point in the reference space
            % xInt: location of GP in the physical space
            % elem: istance of elements class for Basis Function evaluation
            xiMaster = zeros(size(xi,1),2);
            tol = 1e-8;
            itMax = 10;
            nodeM = obj.intMaster.surfaces(elM,:);
            nodeS = obj.intSlave.surfaces(elS,:);
            coord = obj.intMaster.coordinates;
            % loop trough each integration point to project
            for ii = 1:size(xi,1)
                iter = 0;
                w = 0;
                Nslave = elemSlave.computeBasisF(xi(ii,:));
                Nmaster = elemMaster.computeBasisF(xiMaster(ii,:));
                % evaluate normal at integration point
                rhs1 = coord(nodeM,:)'*Nmaster';
                nS = n(nodeS,:)'*Nslave';
                rhs = rhs1 - w*nS - xInt(ii,:)';

                % project Gauss Point
                while (norm(rhs,2) > tol) && (iter < itMax)
                   % compute Jacobian of the projection
                   iter = iter + 1;
                   % get local derivatives of basis functions
                   dN = elemMaster.computeDerBasisF(xiMaster(ii,:));
                   J1 = dN*coord(nodeM,:);
                   J = [J1', -nS];
                   ds = J\(-rhs);
                   xiMaster(ii,:) = xiMaster(ii,:) + ds(1:2)';
                   w = w + ds(3);
                   Nmaster = elemMaster.computeBasisF(xiMaster(ii,:));
                   % evaluate normal at integration point
                   rhs1 = coord(nodeM,:)'*Nmaster';
                   rhs = rhs1 - w*nS - xInt(ii,:)';
                end

                if iter >= itMax
                   xiMaster(ii) = nan;
                   % mark gauss point that did not converge
                end
            end
        end
        %

        function [n_n,varargout] = computeNodalNormal(obj,elem)
           % Return vector n of weighed normals
           n_n = zeros(length(obj.intSlave.coordinates),3);
           area = elem.findAreaAndCentroid(1:obj.intSlave.nSurfaces); % area of each cell
           elem_normal = elem.computeNormal(1:obj.intSlave.nSurfaces);
           obj.elemNormal = elem_normal;
           topol = obj.intSlave.surfaces;
           for i = 1:length(obj.intSlave.coordinates)
              % get elements sharing node i
              elems = find(any(ismember(topol,i),2));
              n_n(i,:) = (area(elems)'*elem_normal(elems,:))/sum(area(elems));
           end
           n_n = n_n./vecnorm(n_n,2,2);
           if nargout > 1
              areaNod = elem.computeAreaNod(obj.intSlave);
              n_a = n_n.*areaNod;
              varargout{1} = n_a;
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
            fM = f(obj.masterCoord(:,1),obj.masterCoord(:,2)); % analytical function computed on master mesh
            fS = f(obj.slaveCoord(:,1),obj.slaveCoord(:,2)); % analytical function computed on master mesh
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
            r = obj.computeRBFradius(ptsInt);
            d = sqrt((ptsInt(:,1) - ptsInt(:,1)').^2 + (ptsInt(:,2) - ptsInt(:,2)').^2 + (ptsInt(:,3) - ptsInt(:,3)').^2);
            fiMM = obj.rbfInterp(d,r,type);
        end

        function [Nm,id] = getMasterBasis(obj,elemId,ptsGauss)
           % elemId: iD of the master elements in the interface mesh
           % ptsGauss: real coordinates of gauss points
           % Nm: master basis function matrix
           % id -> logical array with true = active Gauss points
           tol = 1e-4;
           ptsInt = obj.ptsRBF(:,repNum(3,elemId));
           [fiNM,id1] = obj.computeRBFfiNM(ptsInt,ptsGauss,'gauss');
           switch obj.degree
              case 1
                 Nm = (fiNM*obj.wF(:,repNum(obj.nNmaster,elemId)))./(fiNM*obj.w1(:,elemId));
                 Nsupp = Nm(:,[1 2 3]);
              case 2
                 Ntmp = (fiNM*obj.wF(:,repNum(mortar.nNmaster+2,elemId)))./(fiNM*obj.w1(:,elemId));
                 Nm = Ntmp(:,1:mortar.nNmaster);
                 Nsupp = Nm(:,[end-1 end]);
           end
           % automatically detect supports computing interpolant
           id = all([Nsupp >= 0-tol id1],2);
        end

        function [fiNM,id] = computeRBFfiNM(obj,ptsInt, ptsGauss, type)
            % id: id of points that has a distance < r with at least one
            % master points
            d = sqrt((ptsGauss(:,1) - ptsInt(:,1)').^2 + (ptsGauss(:,2) - ptsInt(:,2)').^2 + (ptsGauss(:,3) - ptsInt(:,3)').^2);
            r = sqrt((max(ptsInt(:,1)) - min(ptsInt(:,1)))^2 + (max(ptsInt(:,2)) - min(ptsInt(:,2)))^2 + (max(ptsInt(:,3)) - min(ptsInt(:,3)))^2);
            id = ~all(d>=r,2);
            fiNM = obj.rbfInterp(d,r,type);
        end

        function [wF,w1] = getSuppWeight(obj)
            wF = zeros(3,obj.nElMaster*4); % weights for support interesection
            w1 = zeros(3,obj.nElMaster*2); % weight for rescaling of support interesection
            n1 = [1 2 3];
            n2 = [1 4 3];
            b = [0 1; 0 0; 1 0];
            % loop trough master elements and interpolate Basis function
            for i = 1:obj.nElMaster
                ptsInt = obj.masterCoord(obj.masterTopol(i,:),:); % get real coordinates of master nodes
                fiMM = obj.computeRBFfiMM(ptsInt,'gauss');
                % solve local system to get weight of interpolant
                wF(:,[4*i-3 4*i-2 4*i-1 4*i]) = [fiMM(n1,n1)\b fiMM(n2,n2)\b]; 
                w1(:,[2*i-1 2*i]) = [fiMM(n1,n1)\ones(3,1) fiMM(n2,n2)\ones(3,1)]; 
            end
        end
        % 
        % function id = detectSupports(obj,iM,wFMat,w1Mat,N) 
        %     % detect slave points in ptsG lying in the support of master element iM
        %     switch obj.degree
        %         case 1 % simply use basis F already computed
        %             id = all([N(:,[1 3]) >= 0, N(:,[1 3]) <= 1],2);
        %         case 2 % use auxiliary basis functions on pair of opposite nodes
        %             wF = wFMat(:,iM*obj.nNelem)
        %             NSupp = [(fiNM(:,n1)*wF(:,[1 2]))./(fiNM(:,n1)*w1(:,1)),...
        %     %     (fiNM(:,n2)*wF(:,[3 4]))./(fiNM(:,n2)*w1(:,2))];
        % 
        %     end

            % ptsMaster = obj.masterCoord(obj.masterTopol(iM,:),:);
            % n1 = [1 2 3];
            % n2 = [1 4 3];
            % tol = -1e-2;
            % fiNM = obj.computeRBFfiNM(ptsMaster,ptsG,'gauss');
            % NSupp = [(fiNM(:,n1)*wF(:,[1 2]))./(fiNM(:,n1)*w1(:,1)),...
            %     (fiNM(:,n2)*wF(:,[3 4]))./(fiNM(:,n2)*w1(:,2))];
            % %Nsupp = (fiNM*wF)./(fiNM*w1); 
            % id1 = ~any(NSupp<tol,2);
            % id2 = all([N>=0, N<=1],2);
            % id = all([id1,id2],2);
            %id = all([Nsupp >= 0, Nsupp <= 1],2);
        % end

        

        function getWeights(obj,interface,nInt,elem,type)
           switch obj.masterCellType
              case 10
                 numPts = sum(1:nInt);
              case 12
                 numPts = (nInt)^2;
           end
           switch interface
              case 'master'
                 switch obj.degree
                    case 1
                       weighF = zeros(numPts,obj.nElMaster*obj.nNmaster);
                       weigh1 = zeros(numPts,obj.nElMaster);
                       pts = zeros(numPts,obj.nElMaster*3);
                       for i = 1:obj.nElMaster
                          [f, ptsInt] = computeMortarBasisF(obj,i, nInt, obj.masterTopol, obj.masterCoord, elem);
                          fiMM = obj.computeRBFfiMM(ptsInt,type);
                          % solve local system to get weight of interpolant
                          weighF(:,repNum(obj.nNmaster,i)) = fiMM\f;
                          weigh1(:,i) = fiMM\ones(size(ptsInt,1),1);
                          pts(:,repNum(3,i)) = ptsInt;
                       end
                    case 2 % second order interpolation (8 nodes + 2 aux)
                       weighF = zeros(numPts,obj.nElMaster*(obj.nNmaster+2));
                       weigh1 = zeros(numPts,obj.nElMaster);
                       pts = zeros(numPts,obj.nElMaster*3);
                       for i = 1:obj.nElMaster
                          [f, ptsInt,fL] = computeMortarBasisF(obj,i, nInt, obj.masterTopol, obj.masterCoord, elem);
                          fiMM = obj.computeRBFfiMM(ptsInt,type);
                          % solve local system to get weight of interpolant
                          weighF(:,repNum(10,i)) = fiMM\[f,fL]; % 8 basis functions + 2 auxiliary for support detection
                          weigh1(:,i) = fiMM\ones(size(ptsInt,1),1);
                          pts(:,repNum(3,i)) = ptsInt;
                       end
                 end
              case 'slave' % to get identity in the conforming case, one should interpolate also the slave basis functions
                 weighF = zeros(numPts,obj.nElSlave*4);
                 weigh1 = zeros(numPts,obj.nElSlave);
                 pts = zeros(numPts,obj.nElSlave*3);
                 for i = 1:obj.nElSlave
                    [f, ptsInt] = computeBasisF2D(i, nInt, obj.slaveTopol, obj.slaveCoord, elem);
                    fiMM = obj.computeRBFfiMM(ptsInt,type);
                    % solve local system to get weight of interpolant
                    weighF(:,[4*i-3 4*i-2 4*i-1 4*i]) = fiMM\f;
                    weigh1(:,i) = fiMM\ones(size(ptsInt,1),1);
                    pts(:,[3*i-2 3*i-1 3*i]) = ptsInt;
                 end
           end
           obj.wF = weighF;
           obj.w1 = weigh1;
           obj.ptsRBF = pts;
           % loop trough master elements and interpolate Basis function
        end

        function elem = getElem(obj,gauss,side)
           % get istance of element class based on cell type
           switch side
              case 'master'
                 elem_tmp = Elements(obj.intMaster,gauss);
                 if obj.intMaster.surfaceVTKType(1) == 5
                    elem = elem_tmp.tri;
                 elseif obj.intMaster.surfaceVTKType(1) == 9
                    elem = elem_tmp.quad;
                 end
              case 'slave'
                 elem_tmp = Elements(obj.intSlave,gauss);
                 if obj.intSlave.surfaceVTKType(1) == 5
                    elem = elem_tmp.tri;
                 elseif obj.intSlave.surfaceVTKType(1) == 9
                    elem = elem_tmp.quad;
                 end
           end
        end

        function E = getMortarOperator(obj,varargin)
           if nargin > 2
              error('Too many input arguments \n');
           end
           E = obj.D\obj.M;
           if isempty(varargin)
              % directly return mortar operator (scalar interpolation)
              return
           else
              nc = varargin{1};
              idRow = 1:size(E,1);
              idCol = 1:size(E,2);
              Enew = sparse(nc*numel(idRow),3*numel(idCol));
              for k = nc-1:-1:0
                 r = nc*idRow - k;
                 c = nc*idCol - k;
                 Enew(r,c) = E;
              end
              E = Enew;
              %
           end
        end

    end

    methods (Access=private)        
       function [bf,pos,varargout] = computeMortarBasisF(obj,elemID,nInt,topol,coord,elem)
          % evaluate shape function in the real space and return position of
          % integration points in the real space (extended to x,y)
          % already ordered to perform RBF interpolation
          %
          % place interpolation points in a regular grid 
          intPts = getIntPoint(obj,nInt);
          bf = computeBasisF(elem,intPts);
          % get coords of interpolation points in the real space
          pos = bf*coord(topol(elemID,:),:);

          % get basis functions of lower order element for support detection (only
          % for quadrilateral elements)
          if nargout==3
             assert(obj.degree==2,'Incorrect number of outputs for basis functions of degree 1')
             bfL = computeBasisF(elem.quadL,intPts);
             varargout{1} = bfL(:,[1 3]);
          end
       end

       function id = detectSupportEB(obj,xi)
          tol = 1e-6;
          switch obj.masterCellType
             case 10
                id = all([xi >= 0-tol, xi <= 1+tol],2);
             case 12
                id = all([xi >= -1-tol, xi <= 1+tol],2);
          end
       end

       function intPts = getIntPoint(obj,nInt)
          switch obj.masterCellType
             case 10
                % uniform grid in triangle
                intPts = zeros(sum(1:nInt),2);
                p = linspace(0,1,nInt);
                k = nInt;
                c = 0;
                for i = 1:nInt
                   intPts(c+1:c+k,1) = p(1:k)';
                   c = c + k;
                   k = k-1;
                end
                intPts(:,2) = repelem(p,nInt:-1:1);              
             case 12
                intPts = linspace(-1,1, nInt);
                [y, x] = meshgrid(intPts, intPts);
                intPts = [x(:), y(:)];
          end
       end

       function connMat = cartConn(obj,mshMaster,mshSlave)
          % get element connectivity for  special case of
          % cartesian meshes in a plane (the procedure is much cheaper)
          cs = mshSlave.coordinates;
          cm = mshMaster.coordinates;
          nS = mshSlave.nSurfaces;
          nM = mshMaster.nSurfaces;
          [xmin,xmax,ymin,ymax] = deal(min(cs(:,1)),max(cs(:,1)),min(cs(:,2)),max(cs(:,2)));
          n0 = round(max(nM,nS)*(max(nM,nS)/min(nM,nS))); % reasonable size for preallocation
          sX = abs(xmax-xmin);
          sY = abs(ymax-ymin);
          dsX = sX/sqrt(nS);
          dsY = sY/sqrt(nS);
          mVec = zeros(n0,1);
          sVec = zeros(n0,1);
          s = 0;
          for im = 1:mshMaster.nSurfaces
             ctmp = cm(mshMaster.surfaces(im,:),1:2);
             [x1,x2] = deal(min(ctmp(:,1)),max(ctmp(:,1)));
             [y1,y2] = deal(min(ctmp(:,2)),max(ctmp(:,2)));
             i1x = fix((x1-xmin)/dsX)+1;
             i2x = ceil((x2-xmin)/dsX);
             lX = i1x:i2x;
             i1y = fix((y1-ymin)/dsY)+1;
             i2y = ceil((y2-ymin)/dsX);
             lY = i1y:i2y;
             [lX,lY] = meshgrid(lX,lY);
             el = sub2ind([sqrt(nS) sqrt(nS)],lX,lY);
             n = numel(el);
             mVec(s+1:s+n) = repelem(im,n);
             sVec(s+1:s+n) = el(:);
             s = s + n;
          end
          mVec = mVec(1:s);
          sVec = sVec(1:s);
          connMat = sparse(mVec,sVec,true(s,1),nM,nS);
       end

       function setupEdgeTopology(obj)
          % reorder surface topology
          % inspired by face topology for FV in MRST
          bot = obj.intSlave.surfaces(:, [1, 2]);
          top = obj.intSlave.surfaces(:, [3, 4]);
          lft = obj.intSlave.surfaces(:, [4, 1]);
          rgt = obj.intSlave.surfaces(:, [2, 3]);
          % unique matrix containing all existing edges (with repetitions)
          edgeMat = [bot;top;lft;rgt];
          %
          [edgeMat,i] = sort(edgeMat,2);
          i = i(:,1);
          %
          %
          %
          % id is [cellnumber, half-face tag]
          id       = [(1:obj.intSlave.nSurfaces)', repmat(1, [obj.intSlave.nSurfaces, 1]);...
                      (1:obj.intSlave.nSurfaces)', repmat(2, [obj.intSlave.nSurfaces, 1]);...
                      (1:obj.intSlave.nSurfaces)', repmat(3, [obj.intSlave.nSurfaces, 1]);...
                      (1:obj.intSlave.nSurfaces)', repmat(4, [obj.intSlave.nSurfaces, 1])];
          % Sort rows to find pairs of cells sharing a face
          [edgeMat, j] = sortrows(edgeMat);
          
          % encode edge matrix
          [obj.edge2nodes,n] = obj.rle(edgeMat);
          obj.nEdges = numel(n);
          N = repelem(1:obj.nEdges,n);
          obj.edge2cells = accumarray([N',i(j)],id(j,1),[obj.nEdges,2]);
       end

       function [mf2c,sf2c] = buildFace2CellMap(obj)
          % get cell ID containing a face
          mf2c = zeros(obj.nElMaster,1);
          sf2c = zeros(obj.nElSlave,1);
          listCell = (1:obj.mshSlave.nCells)';
          for i = 1:obj.nElSlave
             idMat = ismember(obj.mshSlave.cells(listCell,:),obj.slaveTopol(i,:));
             idMat = sum(idMat,2);
             id = find(idMat==4);
             assert(isscalar(id),'Invalid connectivity between face and cells');
             sf2c(i) = listCell(id);
             listCell(id) = [];
          end
          listCell = (1:obj.mshMaster.nCells)';
          for i = 1:obj.nElMaster
             idMat = ismember(obj.mshMaster.cells(listCell,:),obj.masterTopol(i,:));
             idMat = sum(idMat,2);
             id = find(idMat==4);
             assert(isscalar(id),'Invalid connectivity between face and cells');
             mf2c(i) = listCell(id);
             listCell(id) = [];
          end
       end
    end






    methods(Static)
       function [outMat,count] = rle(inMat)
          % the input matrix has already been sorted by rows
          % Find unique rows and their first occurrences
          [outMat, firstIdx, ~] = unique(inMat, 'rows', 'stable');

          % Compute counts of each unique row
          count = diff([firstIdx; size(inMat,1) + 1]);

          % Restore original order of indices
          %indices = sortIdx(firstIdx);
       end
       
        function Nslave_dual = computeDualBasisF(Nslave,gpW)
            % compute dual basis function on slave GPs
            % get matrix of BF phi = Nslave*A
            % A = M\D
            M = Nslave'*(Nslave.*gpW');
            D = diag(Nslave'*gpW');
            A = M\D;
            Nslave_dual = Nslave*A;
        end

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

