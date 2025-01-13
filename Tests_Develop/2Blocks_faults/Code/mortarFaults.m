classdef mortarFaults < handle
   properties
      nDom = 2
      nInterf = 1
      meshGlue    % instance of the meshGlue class (contains mortar algorithms)
      % MORTAR MATRICES
      Dg          % mortar slave matrix in global coords
      Mg          % mortar cross matrix in global coords
      Dn          % mortar slave matrix for normal gaps
      Mn          % mortar cross matrix for normal gaps
      Dt          % mortar slave matrix for tangential gaps
      Mt          % mortar cross matrix for tangential gaps
      L           % lagrange multiplier mass matrix (local coordinates)
      nNslave     % number of nodes per slave elements
      nNmaster    % number
      indN        % index of basis functions in displacement basis matrix
      nS          % number of interface master nodes
      nM          % number of interface slave nodes
      coes        % coesion
      phi         % friction angle
      tagMaster   % id of master domain
      tagSlave     % id of slave domain
      idMaster
      idSlave
      totNodMaster % total number of master domain nodes
      totNodSlave  % total number of slave domain nodes
      nDoF
      activeSet   % group nodes based on active set state
      gap         % current gap vector
      gapOld     % gap of previous time step
      g_T         % \Delta_n g_T
      g_N         % n \cdot g
      E           % mortar operator
      R           % global rotation matrix
      tolGap = 1e-8; % tolerance on normal and tangential gaps
      tolNormal = 1e-3 % tolerance on normal traction
      tolTang = 1e-2 % relative tolerance w.r.t tau_lim
   end

   methods
      function obj = mortarFaults(meshGlue,coes,phi)
         obj.coes = coes;
         obj.phi = phi;
         obj.setParams(meshGlue);
      end

      function setParams(obj,mG)
         obj.meshGlue = mG;
         obj.nNmaster = mG.interfaces.mortar.nNmaster;
         obj.nNslave = mG.interfaces.mortar.nNslave;
         assert(numel(mG.model)==2,'Too many input domains');
         % pattern of basis functions transform
         obj.indN = repmat([1 5 9],1,obj.nNslave)+9*repelem(0:obj.nNslave-1,1,3);
         % save model metrics
         obj.tagMaster = obj.meshGlue.interfaces(1).Master;
         obj.tagSlave = obj.meshGlue.interfaces(1).Slave;
         obj.totNodMaster = obj.meshGlue.model(obj.tagMaster).Grid.topology.nNodes;
         obj.totNodSlave = obj.meshGlue.model(obj.tagSlave).Grid.topology.nNodes;
         obj.nS = length(mG.interfaces.mortar.nodesSlave);
         obj.nM = length(mG.interfaces.mortar.nodesMaster);
         obj.idMaster = obj.meshGlue.interfaces(1).masterSet;
         obj.idSlave = obj.meshGlue.interfaces(1).slaveSet;
         % compute global rotation matrix
         obj.R = getGlobalRotationMatrix(obj);
         % initialize stress (not flexible at all)
         % build dof map with initial active set status
         %
         obj.nDoF = 3*obj.totNodMaster+3*obj.totNodSlave+3*obj.nS;
         obj.E = obj.meshGlue.interfaces.mortar.getMortarOperator(3); % mortar operator
      end


      function computeContactMatrices(obj)
         % get and normalize nodal normals
         area_nod = sqrt(sum(obj.meshGlue.interfaces(1).nodeNormal.^2,2)); % nodal area
         n = obj.meshGlue.interfaces(1).nodeNormal./area_nod;
         % compute mortar matrices within one simulation loop use radial
         % basis functions to evaluate matrix M
         mortar = obj.meshGlue.interfaces(1).mortar; % only one interface defined
         % detect elements lying in the boundary of the interface
         tol = 1e-3;
         boundElem = zeros(mortar.nElSlave,1);
         s =0;
         for i = 1:mortar.nElSlave
            nList = mortar.intSlave.surfaces(i,:);
            yMin = 0; yMax = 10; zMin = 0; zMax = 15;
            idY = any([abs(mortar.intSlave.coordinates(nList,2)-yMin)<tol,abs(mortar.intSlave.coordinates(nList,2)-yMax)<tol],2);
            idZ = any([abs(mortar.intSlave.coordinates(nList,3)-zMin)<tol,abs(mortar.intSlave.coordinates(nList,3)-zMax)<tol],2);
            if any(any([idY,idZ],2))
               boundElem(s+1) = i;
               s = s+1;
            end
         end
         boundElem = boundElem(boundElem~=0);
         %
         nGP = obj.meshGlue.interfaces(1).nG; % numb. of gauss points for element based integration
         mult_type = obj.meshGlue.interfaces(1).multType;
         c_ns = 0;  % counter for GP not projected
         %Mdetect = zeros(mortar.nElMaster,mortar.nElSlave);
         % set Gauss class
         gM = Gauss(mortar.masterCellType,3,2); % gauss class for Master element interpolation
         gS = Gauss(mortar.slaveCellType,nGP,2); % gauss class for slave integration
         elemSlave = getElem(mortar,gS,'slave');
         [imVec,jmVec,Mgvec,Mtvec,Mnvec] = deal(zeros(nnz(mortar.elemConnectivity)*mortar.nNmaster^2,1));
         [isVec,jsVec,Dgvec,Dtvec,Dnvec] = deal(zeros(nnz(mortar.elemConnectivity)*mortar.nNmaster^2,1));
         % Perform interpolation on the master side (computing weights
         % and interpolation coordinates)
         % Loop trough slave elements
         cs = 0; % slave matrix entry counter
         cm = 0; % master matrix entry counter
         for j = 1:mortar.nElSlave
            %Compute Slave quantities
            dJWeighed = elemSlave.getDerBasisFAndDet(j,3); % Weighted Jacobian
            %get Gauss Points position in the real space
            ptsGauss = getGPointsLocation(elemSlave,j);
            nSlave = mortar.intSlave.surfaces(j,:);
            n_el = (n(nSlave,:))';
            n_el = n_el(:); % local nodal normal vector
            Rloc = getRotationMatrix(obj,n_el);
            %Rloc = eye(3*obj.nNslave);
            %A = diag(repelem(area_nod(nSlave),3));
            NSlave = getBasisFinGPoints(elemSlave); % Get slave basis functions
            switch mult_type
               case 'standard'
                  NSlaveMult = NSlave; % Slave basis functions
               case 'dual'
                  NSlaveMult = mortar.computeDualBasisF(NSlave,dJWeighed);
            end
            % modify multiplier basis for elements containing boundary
            % nodes 
            % if ismember(j,boundElem)
            %    NSlaveMult = ones(size(NSlaveMult))/size(NSlaveMult,2); 
            %    % constant unit function
            % end
            master_elems = find(mortar.elemConnectivity(:,j));
            for jm = master_elems'
               nMaster = mortar.intMaster.surfaces(jm,:);
               [NMaster,id] = mortar.getMasterBasis(jm,ptsGauss); % compute interpolated master basis function \Pi(Nm)
               if any(id)
                  % element-based integration
                  % prepare provisional 3D matrices
                  Nm = obj.dispSP(NMaster(id,:));
                  Ns = obj.dispSP(NSlave(id,:));
                  Nmult = obj.dispSP(NSlaveMult(id,:));
                  % global mortar matrices (global frame)
                  Dgtmp = pagemtimes(Nmult,'transpose',Ns,'none');
                  Dgtmp = Dgtmp.*reshape(dJWeighed(id),1,1,[]);
                  Dgloc = Rloc'*sum(Dgtmp,3);
                  Mgtmp = pagemtimes(Nmult,'transpose',Nm,'none');
                  Mgtmp = Mgtmp.*reshape(dJWeighed(id),1,1,[]);
                  Mgloc = Rloc'*sum(Mgtmp,3);
                  % normal mortar matrices (global frame)
                  Nn = pagemtimes(Ns,n_el);
                  Dntmp = pagemtimes(Nmult(1,:,:),'transpose',pagemtimes(Nn,'transpose',Ns,'none'),'none');
                  Dntmp = Dntmp.*reshape(dJWeighed(id),1,1,[]);
                  Dnloc = Rloc'*sum(Dntmp,3);
                  Mntmp = pagemtimes(Nmult(1,:,:),'transpose',pagemtimes(Nn,'transpose',Nm,'none'),'none');
                  Mntmp = Mntmp.*reshape(dJWeighed(id),1,1,[]);
                  Mnloc = Rloc'*sum(Mntmp,3);
                  % tangential mortar matrices (global frame)
                  Nt = eye(3) - pagemtimes(Nn,'none',Nn,'transpose');
                  Dttmp = pagemtimes(Nmult,'transpose',pagemtimes(Nt,Ns),'none');
                  Dttmp = Dttmp.*reshape(dJWeighed(id),1,1,[]);
                  Dtloc = Rloc'*sum(Dttmp,3);
                  Mttmp = pagemtimes(Nmult,'transpose',pagemtimes(Nt,Nm),'none');
                  Mttmp = Mttmp.*reshape(dJWeighed(id),1,1,[]);
                  Mtloc = Rloc'*sum(Mttmp,3);
                  % lagrange multipliers mass matrix L (local frame!)
                  Ltmp = pagemtimes(Nmult,'transpose',Nmult,'none');
                  Ltmp = Ltmp.*reshape(dJWeighed(id),1,1,[]);
                  Lloc = sum(Ltmp,3);
                  % local assembly
                  dof_master = get_dof(nMaster);
                  dof_slave = get_dof(nSlave);
                  [jjM,iiM] = meshgrid(dof_master,dof_slave);
                  [jjS,iiS] = meshgrid(dof_slave,dof_slave);
                  nm = numel(Mgloc);
                  ns = numel(Dgloc);
                  imVec(cm+1:cm+nm) = iiM(:); jmVec(cm+1:cm+nm) = jjM(:);
                  isVec(cs+1:cs+ns) = iiS(:); jsVec(cs+1:cs+ns) = jjS(:);
                  Mgvec(cm+1:cm+nm) = Mgloc(:);
                  Dgvec(cs+1:cs+ns) = Dgloc(:);
                  Mnvec(cm+1:cm+nm) = Mnloc(:);
                  Dnvec(cs+1:cs+ns) = Dnloc(:);
                  Mtvec(cm+1:cm+nm) = Mtloc(:);
                  Dtvec(cs+1:cs+ns) = Dtloc(:);
                  Lvec(cs+1:cs+ns)  = Lloc(:);
                  % sort out Points already projected
                  dJWeighed = dJWeighed(~id);
                  ptsGauss = ptsGauss(~id,:);
                  NSlave = NSlave(~id,:);
                  NSlaveMult = NSlaveMult(~id,:);
                  cs = cs+ns;
                  cm = cm+nm;
               end
            end
            if ~all(id)
               fprintf('GP not sorted for slave elem %i \n',j);
               c_ns = c_ns + 1;
            end
         end
         imVec = imVec(1:cm); jmVec = jmVec(1:cm);
         isVec = isVec(1:cs); jsVec = jsVec(1:cs);
         Mgvec = Mgvec(1:cm); Mnvec = Mnvec(1:cm); Mtvec = Mtvec(1:cm);
         Dgvec = Dgvec(1:cm); Dnvec = Dnvec(1:cm); Dtvec = Dtvec(1:cm);
         Lvec = Lvec(1:cs);
         obj.Mg = sparse(imVec,jmVec,Mgvec,3*obj.nS,3*obj.nM);
         obj.Mn = sparse(imVec,jmVec,Mnvec,3*obj.nS,3*obj.nM);
         obj.Mt = sparse(imVec,jmVec,Mtvec,3*obj.nS,3*obj.nM);
         obj.Dg = sparse(isVec,jsVec,Dgvec,3*obj.nS,3*obj.nS);
         obj.Dn = sparse(isVec,jsVec,Dnvec,3*obj.nS,3*obj.nS);
         obj.Dt = sparse(isVec,jsVec,Dtvec,3*obj.nS,3*obj.nS);
         obj.L = sparse(isVec,jsVec,Lvec,3*obj.nS,3*obj.nS);
         % eliminate small numerical quantities in obj.Dg
         %obj.Dg(abs(obj.Dg)<1e-12) = 0;
      end

      function R = getRotationMatrix(obj,n_el)
         % n_el: local nodal normal array
         nNode = round(numel(n_el)/3); 
         R = zeros(3*nNode,3*nNode);
         for i = 0:nNode-1
            n = n_el(3*i+1:3*i+3);
            [Rn,~] = qr(n);
            v = Rn'*n;
            if v(1)>0
               Rn(:,1) = -Rn(:,1);
            end
            R(3*i+1:3*i+3,3*i+1:3*i+3) = Rn;
         end
      end

      function R = getNodeRotationMatrix(obj,k)
         % n_el: local nodal normal array
         n = obj.meshGlue.interfaces.nodeNormal(k,:);
         n = (n/norm(n,2))';
         [R,~] = qr(n);
         v = R'*n;
         if v(1)>0
            R(:,1) = -R(:,1);
         end
      end
      

      function Rloc = getGlobalRotationMatrix(obj)
         area_nod = sqrt(sum(obj.meshGlue.interfaces(1).nodeNormal.^2,2));
         n_el = obj.meshGlue.interfaces(1).nodeNormal./area_nod; % normalized nodal area
         % n_el: local nodal normal array
         nEntry = 9;
         Rvec = zeros(nEntry*obj.nS,1);
         iivec = zeros(nEntry*obj.nS,1);
         jjvec = zeros(nEntry*obj.nS,1);
         l1 = 0;
         for i = 0:obj.nS-1
            n = (n_el(i+1,:))';
            [Rloc,~] = qr(n);
            v = Rloc'*n;
            if v(1)>0
               Rloc(:,1) = -Rloc(:,1);
            end
            [iiLoc,jjLoc] = meshgrid(3*i+1:3*i+3,3*i+1:3*i+3);
            iivec(l1+1:l1+nEntry) = iiLoc(:);
            jjvec(l1+1:l1+nEntry) = jjLoc(:);
            Rvec(l1+1:l1+nEntry) = Rloc(:);
            l1 = l1+nEntry;
         end
         Rloc = sparse(iivec,jjvec,Rvec,3*obj.nS,3*obj.nS);
      end

      function Nout = dispSP(obj,Nin)
         % reshape basis functions to obtain displacement shape function
         % input: nG x nN matrix
         % output: 3 x nN x nG
         s2 = 3*size(Nin,2);
         s3 = size(Nin,1);
         N = zeros(3*s2,s3);
         Nrep = repelem(Nin',3,1);
         N(obj.indN(1:s2),:) = Nrep;
         Nout = reshape(N,[3,s2,s3]); % reshaped 3D matrix
      end

      function dofs = getContactDofs(obj,tag,varargin)
         if ~isempty(varargin)
            dofIn = varargin{1};
            dofIn = reshape(dofIn,1,[]);
         else
            dofIn = [];
         end
         % dofIn must be a row vector
         % map nodal entries to degrees of freedom in the global linear
         % system. 
         % tag: 'master','slave','lag'
         switch tag
            case 'master'
               if isempty(dofIn)
                  dofIn = 1:obj.totNodMaster;
               end
               dofs = get_dof(dofIn);
            case 'slave'
               if isempty(dofIn)
                  dofIn = 1:obj.totNodSlave;
               end
               dofs = 3*obj.totNodMaster+get_dof(dofIn);
            case 'lag'
               if isempty(dofIn)
                  dofs = [];
                  return
               end
               dofs = 3*(obj.totNodMaster+obj.totNodSlave)+get_dof(dofIn);
            otherwise
               error(['Invalide tag string for getContactDofs method \n:' ...
                  'valid inputs are: master, slave, lag']);
         end
      end

      function tLim = computeLimitTraction(obj,activeSet,mult)
         % TO DO: add check on magnitud gap
         % get limit traction array in local coordinates
         tLim = zeros(3*numel(activeSet.curr.slip),1);
         l = 1;
         for i = activeSet.curr.slip'
            dofSlip = get_dof(i);
            gapNorm = norm(obj.g_T(dofSlip),2);
            if gapNorm > obj.tolGap
               % principle of maximum plastic dissipation
               tL = repelem(tau_max(obj,i,mult),3,1).*obj.g_T(dofSlip)/gapNorm;
               tL = obj.getNodeRotationMatrix(obj.idSlave(i))'*tL; % traction in local coords
               tL(1) = 0; 
            else
               %keep direction of previous lagrange multiplier
               lambda = mult(get_dof(i));
               lambda(1) = 0;
               tL = repelem(tau_max(obj,i,mult),3,1).*(lambda/norm(lambda));
               %this is already in local components
            end
            tLim(3*l-2:3*l) = tL;
            l = l+1;
         end
      end

      function tau = tau_max(obj,nList,multipliers)
         % get tau max for input list of nodes
         % use the normal component of contact traction
         tau = obj.coes - tan(deg2rad(obj.phi))*multipliers(3*nList-2);
      end

      function computeNodalGap(obj,state,dofMap)
         % compute gap in global coordinates
         % compute time difference of the tangential component of nodal gap
         % \Delta_n g_T (in global coordinates)
         % Use mortar operator E to map master nodes to slave side
         usCurr = state(obj.tagSlave).dispCurr(dofMap.nodSlave);
         umCurr = state(obj.tagMaster).dispCurr(dofMap.nodMaster);
         obj.gap = usCurr - obj.E*umCurr;
         % compute weighted nodal gap
         gapW = obj.Dg*usCurr - obj.Mg*umCurr;
         if isempty(obj.gapOld)
            obj.gapOld = obj.gap;
         end
         obj.g_T = obj.gap - obj.gapOld; % global nodal gap
         area_nod = sqrt(sum(obj.meshGlue.interfaces.nodeNormal.^2,2)); % nodal area
         n = obj.meshGlue.interfaces.nodeNormal./area_nod; % unit nodal normal
         % compute tangential component
         obj.g_N = zeros(obj.nS,1);
         l1 = 0;
         for i = 1:obj.nS
            % get node normal
            n_i = n(i,:);
            T = eye(3) - n_i'*n_i; % tangential projection matrix
            obj.g_T(l1+1:l1+3) = T*obj.g_T(l1+1:l1+3);
            obj.g_N(i) = n_i*obj.gap(l1+1:l1+3);
            l1 = l1+3;
         end
      end

      function Tmat = computeDtDgt(obj,activeSet,mult)
         % compute derivative of tangential traction w.r.t tangential gap
         % quantities are in global coordinates
         iVec = zeros(9*numel(activeSet.curr.slip),1);
         jVec = zeros(9*numel(activeSet.curr.slip),1);
         Tvec = zeros(9*numel(activeSet.curr.slip),1);
         l = 0;
         s = 1;
         nSlip = numel(activeSet.curr.slip);
         for i = activeSet.curr.slip' 
            dofLoc = get_dof(s);
            dofMult = get_dof(i);
            tLim = tau_max(obj,i,mult);
            DgT = -obj.getNodeRotationMatrix(obj.idSlave(i))*obj.g_T(dofMult);
            normgT = norm(DgT);
            if normgT > obj.tolGap
               Tloc = tLim*(normgT^2*eye(3)-DgT*DgT')/normgT^3; % 3x3 local mat in global coords
               [ii,jj] = meshgrid(dofLoc,dofLoc);
               iVec(l+1:l+9) = ii(:);
               jVec(l+1:l+9) = jj(:);
               Tvec(l+1:l+9) = Tloc(:);
               l = l+9;
            end
            s = s+1;
         end
        iVec = iVec(iVec~=0); 
        jVec = jVec(jVec~=0);
        if isempty(iVec) && isempty(jVec)
           Tmat = sparse(3*nSlip,3*nSlip);
        else
           Tmat = sparse(iVec,jVec,Tvec,3*nSlip,3*nSlip);
        end
      end

      function [TD,TM,N] = computeConsistencyMatrices(obj,activeSet,mult,state,dofMap)
         % compute the non linear consistency matrices by directly evaluating
         % the gap at the gauss points of the slave side
         usCurr = state(obj.tagSlave).dispCurr(dofMap.nodSlave);
         umCurr = state(obj.tagMaster).dispCurr(dofMap.nodMaster);
         area_nod = sqrt(sum(obj.meshGlue.interfaces(1).nodeNormal.^2,2)); % nodal area
         nvec = obj.meshGlue.interfaces(1).nodeNormal./area_nod;
         mortar = obj.meshGlue.interfaces(1).mortar; % only one interface defined
         nGP = obj.meshGlue.interfaces(1).nG; % numb. of gauss points for element based integration (use this for all matrices)
         mult_type = obj.meshGlue.interfaces(1).multType;
         c_ns = 0;  % counter for GP not projected
         %Mdetect = zeros(mortar.nElMaster,mortar.nElSlave);
         % set Gauss class
         gS = Gauss(mortar.slaveCellType,nGP,2); % gauss class for slave integration
         elemSlave = getElem(mortar,gS,'slave');
         [imVec,jmVec,Mtvec] = deal(zeros(nnz(mortar.elemConnectivity)*mortar.nNmaster^2,1));
         [isVec,jsVec,Dtvec] = deal(zeros(nnz(mortar.elemConnectivity)*mortar.nNslave^2,1));
         [inVec,jnVec,Nvec] = deal(zeros(mortar.nElSlave*mortar.nNslave^2,1));
         cn = 0; % counter for normal matrix
         cm = 0;
         cs = 0;
         for j = 1:mortar.nElSlave
            %Compute Slave quantities
            nSlave = mortar.intSlave.surfaces(j,:); % slave node id
            activeNodes = ismember(nSlave,activeSet.curr.slip);
            if ~any(activeNodes)
               continue
            end
            if norm(obj.g_T(get_dof(nSlave))) < obj.tolGap
               % gap too small to compute consistency matrices
               continue
            end
            dJWeighed = elemSlave.getDerBasisFAndDet(j,3); % Weighted Jacobian
            %get Gauss Points position in the real space
            ptsGauss = getGPointsLocation(elemSlave,j);
            n_el = (nvec(nSlave,:))';
            n_el = n_el(:); % local nodal normal vector
            Rloc = getRotationMatrix(obj,n_el);
            Rloc = Rloc(get_dof(find(activeNodes)),get_dof(find(activeNodes)));
            NSlave = getBasisFinGPoints(elemSlave); % Get slave basis functions
            switch mult_type
               case 'standard'
                  NSlaveMult = NSlave; % Slave basis functions
               case 'dual'
                  NSlaveMult = mortar.computeDualBasisF(NSlave,dJWeighed);
            end
            % compute slave only matrices
            Ns = obj.dispSP(NSlave);
            Nmult = obj.dispSP(NSlaveMult(:,activeNodes));
            gapGP = computeGapInGP(obj,mortar,Ns,usCurr,umCurr,j,ptsGauss);
            dtdgt = computeDerTracGap(obj,Ns,Nmult,nSlave(activeNodes),mult,); % 3D matrix with NL derivatives of tractions
            dtdtn = computeDerTracTn(obj,Ns,nSlave(activeNodes),Rloc); % 3D matrix with NL derivatives of tractions
            Ntmp = pagemtimes(Nmult([2 3],:,:),'transpose',pagemtimes(dtdtn,Nmult(1,:,:)),'none');
            Ntmp =  Ntmp.*reshape(dJWeighed,1,1,[]);
            Nloc = sum(Ntmp,3);
            dof_slave = get_dof(nSlave);
            dof_mult = get_dof(nSlave(activeNodes));
            [jjS,iiS] = meshgrid(dof_slave,dof_mult);
            [jjMult,iiMult] = meshgrid(dof_mult,dof_mult);
            nN = numel(Nloc);
            Nvec(cn+1:cn+nN) = Nloc(:);
            jnVec(cn+1:cn+nN) = jjMult(:); inVec(cn+1:cn+nN) = iiMult(:);
            cn = cn + nN;
            % compute cross-grid matrices
            master_elems = find(mortar.elemConnectivity(:,j)); 
            for jm = master_elems'
               nMaster = mortar.intMaster.surfaces(jm,:);
               [NMaster,id] = mortar.getMasterBasis(jm,ptsGauss); % compute interpolated master basis function \Pi(Nm)
               if any(id)
                  % element-based integration
                  % prepare provisional 3D matrices
                  Nm = obj.dispSP(NMaster(id,:));
                  Ns = obj.dispSP(NSlave(id,:));
                  Nmult = obj.dispSP(NSlaveMult(id,activeNodes));
                  dtdgtTmp = dtdgt(:,:,id);
                  % normal mortar matrices (global frame)
                  Nn = pagemtimes(Ns,n_el);
                  % tangential mortar matrices (global frame)
                  Nt = eye(3) - pagemtimes(Nn,'none',Nn,'transpose');
                  tangNs = pagemtimes(Nt,Ns);
                  tangNm = pagemtimes(Nt,Nm);
                  Dttmp = pagemtimes(Nmult,'transpose',pagemtimes(dtdgtTmp,tangNs),'none');
                  Dttmp = Dttmp.*reshape(dJWeighed(id),1,1,[]);
                  Dtloc = Rloc'*sum(Dttmp,3);
                  Mttmp = pagemtimes(Nmult,'transpose',pagemtimes(dtdgtTmp,tangNm),'none');
                  Mttmp = Mttmp.*reshape(dJWeighed(id),1,1,[]);
                  Mtloc = Rloc'*sum(Mttmp,3);
                  % local assembly
                  dof_master = get_dof(nMaster);
                  [jjM,iiM] = meshgrid(dof_master,dof_mult);
                  nm = numel(Mtloc);
                  ns = numel(Dtloc);
                  imVec(cm+1:cm+nm) = iiM(:); jmVec(cm+1:cm+nm) = jjM(:);
                  isVec(cs+1:cs+ns) = iiS(:); jsVec(cs+1:cs+ns) = jjS(:);
                  Mtvec(cm+1:cm+nm) = Mtloc(:);
                  Dtvec(cs+1:cs+ns) = Dtloc(:);
                  % sort out Points already projected
                  dJWeighed = dJWeighed(~id);
                  ptsGauss = ptsGauss(~id,:);
                  NSlave = NSlave(~id,:);
                  NSlaveMult = NSlaveMult(~id,:);
                  dtdgt = dtdgt(:,:,~id);
                  cs = cs+ns;
                  cm = cm+nm;
               end
            end
            if ~all(id)
               fprintf('GP not sorted for slave elem %i \n',j);
               c_ns = c_ns + 1;
            end
         end

         if cn == 0
            TM = sparse(3*obj.nS,3*obj.nM);
            TD = sparse(3*obj.nS,3*obj.nS);
            N = sparse(3*obj.nS,3*obj.nS);
         else
            imVec = imVec(1:cm); jmVec = jmVec(1:cm);
            isVec = isVec(1:cs); jsVec = jsVec(1:cs);
            inVec = inVec(1:cn); jnVec = jnVec(1:cn);
            Mtvec = Mtvec(1:cm); Dtvec = Dtvec(1:cs); Nvec = Nvec(1:cn);
            TM = sparse(imVec,jmVec,Mtvec,3*obj.nS,3*obj.nM);
            TD = sparse(isVec,jsVec,Dtvec,3*obj.nS,3*obj.nS);
            N = sparse(inVec,jnVec,Nvec,3*obj.nS,3*obj.nS);
         end
      end

      function Nmat = computeDtDtn(obj,activeSet,mult)
         % compute derivative of tangential traction w.r.t tangential gap
         % quantities are in local coordinates
         iVec = zeros(2*numel(activeSet.curr.slip),1);
         jVec = zeros(2*numel(activeSet.curr.slip),1);
         Nvec = zeros(2*numel(activeSet.curr.slip),1);
         l = 0;
         s = 1;
         nSlip = numel(activeSet.curr.slip);
         for i = activeSet.curr.slip'
            dofMult = get_dof(i);
            dofLoc = get_dof(s);
            %tLim = tau_max(obj,i,mult);
            DgT = -obj.getNodeRotationMatrix(i)'*obj.g_T(dofMult); % local tangential gap
            DgT = DgT([2 3]); % if evrything is ok, the first component is actually 0
            normgT = norm(DgT);
            %lambda_t = mult(dofMult);
            %lambda_t = lambda_t([2 3]);
            %normLambda = norm(lambda_t);
            if normgT > obj.tolGap
               Nloc = -tan(deg2rad(obj.phi))*DgT/normgT; % 3x3 local mat in global coords
               %Nloc = tLim*(normgT^2*eye(3)-DgT*DgT')/normgT^3;
               iVec(l+1:l+2) = [dofLoc(2);dofLoc(3)];
               jVec(l+1:l+2) = [dofLoc(1);dofLoc(1)];
               Nvec(l+1:l+2) = Nloc(:);
               l = l+2;
            end
            s = s+1;
         end
         iVec = iVec(iVec~=0);
         jVec = jVec(jVec~=0);
         if isempty(iVec) && isempty(jVec)
            Nmat = sparse(3*nSlip,3*nSlip);
         else
            Nmat = sparse(iVec,jVec,Nvec,3*nSlip,3*nSlip);
         end
      end

      function rhs = computeRhsLimitTraction(obj,mult,activeSet)
         % compute rhs related to limit tracion: <mu,t*>_slip
         mortar = obj.meshGlue.interfaces(1).mortar; % only one interface defined
         nGP = obj.meshGlue.interfaces(1).nG; % numb. of gauss points for element based integration (use this for all matrices)
         mult_type = obj.meshGlue.interfaces(1).multType;
         %Mdetect = zeros(mortar.nElMaster,mortar.nElSlave);
         area_nod = sqrt(sum(obj.meshGlue.interfaces(1).nodeNormal.^2,2)); % nodal area
         nvec = obj.meshGlue.interfaces(1).nodeNormal./area_nod; % normalized nodal normals
         % set Gauss class
         gS = Gauss(mortar.slaveCellType,4,2); % gauss class for slave integration
         elemSlave = getElem(mortar,gS,'slave');
         rhs = zeros(length(mult),1);
         % compute nodal limit traction
         tLim = obj.computeLimitTraction(activeSet,mult);
         for j = 1:mortar.nElSlave
            %Compute Slave quantities
            nSlave = mortar.intSlave.surfaces(j,:); % slave node id
            activeNodes = ismember(nSlave,activeSet.curr.slip);
            if ~any(activeNodes)
               continue
            end
            dofSlave = get_dof(nSlave(activeNodes));
            [id,pos] = ismember(activeSet.curr.slip,nSlave(activeNodes));
            [~,p] = sort(pos(id));
            idN = find(id);
            idN = idN(p);
            t = tLim(get_dof(idN));
            %sn = mult(dofSlave(1:3:end));
            tauLimit = tau_max(obj,nSlave(activeNodes),mult);
            dJWeighed = elemSlave.getDerBasisFAndDet(j,3); % Weighted Jacobian
            %get Gauss Points position in the real space
            NSlave = getBasisFinGPoints(elemSlave); % Get slave basis functions
            switch mult_type
               case 'standard'
                  NSlaveMult = NSlave; % Slave basis functions
               case 'dual'
                  NSlaveMult = mortar.computeDualBasisF(NSlave,dJWeighed);
            end
            % prepare 3D matrices for Gauss integration
            Ns = obj.dispSP(NSlave(:,activeNodes));
            Nmult = obj.dispSP(NSlaveMult(:,activeNodes));       
            if norm(obj.g_T(get_dof(nSlave))) < obj.tolGap
               tSlip = pagemtimes(Nmult,t);
            else
               n_el = (nvec(nSlave(activeNodes),:))';
               n_el = n_el(:); % local nodal normal vector
               Rloc = getRotationMatrix(obj,n_el);
               tauLim = pagemtimes(Nmult(1,1:3:end,:),tauLimit);
               gt = Rloc'*obj.g_T(get_dof(nSlave(activeNodes)));
               gtGP = pagemtimes(Ns,gt);    % interpolated gap
               norm_gt = pagenorm(gtGP);
               tSlip = pagemtimes(tauLim,pagemtimes(gtGP,1/norm_gt)); % 3D limit stress - should i use weighed gap instead?
            end
            rhsTmp = pagemtimes(Nmult(2:3,:,:),'transpose',tSlip(2:3,:,:),'none');
            rhsTmp = rhsTmp.*reshape(dJWeighed,1,1,[]);
            rhsLoc = sum(rhsTmp,3);
            rhs(dofSlave) = rhs(dofSlave)+rhsLoc;
         end
      end


      function rhs = computeRhsLimitTraction2(obj,mult,state,activeSet,dofMap)
         % compute rhs related to limit tracion: <mu,t*>_slip
         % gap is comptuted on each gauss point using exact Pi operator
         usCurr = state(obj.tagSlave).dispCurr(dofMap.nodSlave);
         umCurr = state(obj.tagMaster).dispCurr(dofMap.nodMaster);
         mortar = obj.meshGlue.interfaces(1).mortar; % only one interface defined
         nGP = obj.meshGlue.interfaces(1).nG; % numb. of gauss points for element based integration (use this for all matrices)
         mult_type = obj.meshGlue.interfaces(1).multType;
         %Mdetect = zeros(mortar.nElMaster,mortar.nElSlave);
         area_nod = sqrt(sum(obj.meshGlue.interfaces(1).nodeNormal.^2,2)); % nodal area
         nvec = obj.meshGlue.interfaces(1).nodeNormal./area_nod; % normalized nodal normals
         % set Gauss class
         gS = Gauss(mortar.slaveCellType,4,2); % gauss class for slave integration
         elemSlave = getElem(mortar,gS,'slave');
         rhs = zeros(length(mult),1);
         tLim = obj.computeLimitTraction(activeSet,mult);
         for j = 1:mortar.nElSlave
            %Compute Slave quantities
            nSlave = mortar.intSlave.surfaces(j,:); % slave node id
            activeNodes = ismember(nSlave,activeSet.curr.slip);
            if ~any(activeNodes)
               continue
            end
            dofSlave = get_dof(nSlave(activeNodes));
            [id,pos] = ismember(activeSet.curr.slip,nSlave(activeNodes));
            [~,p] = sort(pos(id));
            idN = find(id);
            idN = idN(p);
            t = tLim(get_dof(idN));
            %sn = mult(dofSlave(1:3:end));
            tauLimit = tau_max(obj,nSlave(activeNodes),mult);
            dJWeighed = elemSlave.getDerBasisFAndDet(j,3); % Weighted Jacobian
            %get Gauss Points position in the real space
            ptsGauss = getGPointsLocation(elemSlave,j);
            NSlave = getBasisFinGPoints(elemSlave); % Get slave basis functions
            switch mult_type
               case 'standard'
                  NSlaveMult = NSlave; % Slave basis functions
               case 'dual'
                  NSlaveMult = mortar.computeDualBasisF(NSlave,dJWeighed);
            end
            n_active = (nvec(nSlave(activeNodes),:))';
            n_active = n_active(:); % local nodal normal vector
            Rloc = getRotationMatrix(obj,n_active);
            % prepare 3D matrices for Gauss integration    
            if norm(obj.g_T(get_dof(nSlave))) < obj.tolGap
               Nmult = obj.dispSP(NSlaveMult(:,activeNodes));
               tSlip = pagemtimes(Nmult,t);
            else
               Ns = obj.dispSP(NSlave);
               Nmult = obj.dispSP(NSlaveMult(:,activeNodes));
               n_el = (nvec(nSlave,:))';
               n_el = n_el(:); % local nodal normal vector
               tauLim = pagemtimes(Nmult(1,1:3:end,:),tauLimit);
               % interpolated slave displacement on each gp
               % us = pagemtimes(Ns,usCurr(get_dof(nSlave)));
               % % compute the tangential gap in each gauss point using
               % % mortar projection
               % master_elems = find(mortar.elemConnectivity(:,j));
               % um = zeros(3,1,size(us,3));
               % idM = (1:size(us,3))';
               % for jm = master_elems'
               %    nMaster = mortar.intMaster.surfaces(jm,:);
               %    [NMaster,id] = mortar.getMasterBasis(jm,ptsGauss);
               %    if any(id)
               %       Nm = obj.dispSP(NMaster(id,:));
               %       um(:,:,idM(id)) = pagemtimes(Nm,umCurr(get_dof(nMaster)));
               %    end
               %    ptsGauss = ptsGauss(~id,:);
               %    % index of non projected gp from previous projection
               %    idM = idM(~id);
               % end
               % gapGP = us - um;
               % ptsGauss = getGPointsLocation(elemSlave,j);
               gapGP = computeGapInGP(obj,mortar,Ns,usCurr,umCurr,j,ptsGauss);
               % project gap onto tangential direction
               % normal mortar matrices (global frame)
               Nn = pagemtimes(Ns,n_el);
               % tangential mortar matrices (global frame)
               Nt = eye(3) - pagemtimes(Nn,'none',Nn,'transpose');
               gtGP = pagemtimes(Nt,gapGP);
               norm_gt = pagenorm(gtGP);
               tSlip = pagemtimes(tauLim,pagemtimes(gtGP,1/norm_gt)); nSlave = mortar.intSlave.surfaces(j,:); % slave node id% 3D limit stress - should i use weighed gap instead?
            end
            rhsTmp = pagemtimes(Nmult(2:3,:,:),'transpose',tSlip(2:3,:,:),'none');
            rhsTmp = rhsTmp.*reshape(dJWeighed,1,1,[]);
            rhsLoc = sum(rhsTmp,3);
            rhs(dofSlave) = rhs(dofSlave)+Rloc'*rhsLoc;
         end
      end


   end

   methods (Access=private)
      function dtdgt = computeDerTracGap(obj,Nmult,nodeId,mult,g)
         % N: 3D slave side matrix
         % result 3D matrix of size (3*nN)x(3*nN)xnG
         % nodeId: node list of input element
         nG = size(Nslave,3);
         % gt = obj.g_T(get_dof(nodeId));
         sn = mult(3*nodeId-2);
         dtdgt = repmat(zeros(3,3),1,1,nG);
         for i = 1:nG
            % get limit tau at gauss point
            sigma_n = Nmult(1,1:3:end,i)*sn;
            tauLim = obj.coes - tan(deg2rad(obj.phi))*sigma_n;
            dtdgt(:,:,i) = tauLim*((eye(3)*norm(g)^2 - g*g')/(norm(g))^3);
         end
      end

      function dtdtn = computeDerTracTn(obj,Nslave,nodeId,Rloc)
         nG = size(Nslave,3);
         gt = obj.g_T(get_dof(nodeId));
         % get rotation matrix
         gt = Rloc'*gt; % get tangential gap in local coordinates
         % Nslave, Nmult: 3D slave side matrix
         % result 3D matrix of size (3*nN)x(3*nN)xnG
         dtdtn = repmat(zeros(2,1),1,1,nG);
         for i = 1:nG
            g = Nslave(:,:,i)*gt; % get tangential gap in gauss point
            g = g([2 3]); % keep only tangential component (the first is zero by construction)
            dtdtn(:,:,i) = -tan(deg2rad(obj.phi))*(g/norm(g));
         end
      end

      function gapGP = computeGapInGP(obj,mortar,Ns,usCurr,umCurr,j,ptsGauss)
         nSlave = mortar.intSlave.surfaces(j,:); % slave node id
         us = pagemtimes(Ns,usCurr(get_dof(nSlave)));
         % compute the tangential gap in each gauss point using
         % mortar projection
         master_elems = find(mortar.elemConnectivity(:,j));
         um = zeros(3,1,size(us,3));
         idM = (1:size(us,3))';
         for jm = master_elems'
            nMaster = mortar.intMaster.surfaces(jm,:);
            [NMaster,id] = mortar.getMasterBasis(jm,ptsGauss);
            if any(id)
               Nm = obj.dispSP(NMaster(id,:));
               um(:,:,idM(id)) = pagemtimes(Nm,umCurr(get_dof(nMaster)));
            end
            ptsGauss = ptsGauss(~id,:);
            % index of non projected gp from previous projection
            idM = idM(~id);
         end
         gapGP = us - um;
      end
   end
end
