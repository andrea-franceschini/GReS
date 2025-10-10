classdef mortarFaultsP0 < handle
   % implementing faults with mortar method.
   % P0 (piecewise constant) multipliers with pressure-jump stabilization
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
      nNmaster    % numberof nodes per master elements
      indN        % index of basis functions in displacement basis matrix
      nS          % number of interface master nodes
      nM          % number of interface slave nodes
      nMult       % number of multipliers
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
      R           % global rotation matrix
      tolGap = 1e-8; % tolerance on normal and tangential gaps
      tolNormal = 1e-4 % tolerance on normal traction
      tolTang = 1e-3 % relative tolerance w.r.t tau_lim
      E
      areaMap
      stabMat
   end

   methods
      function obj = mortarFaultsP0(meshGlue,coes,phi)
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
         obj.nMult = mG.interfaces.mortar.nElSlave;
         obj.idMaster = obj.meshGlue.interfaces(1).masterSet;
         obj.idSlave = obj.meshGlue.interfaces(1).slaveSet;
         % compute global rotation matrix
         obj.R = getGlobalRotationMatrix(obj);
         % initialize stress (not flexible at all)
         % build dof map with initial active set status
         %
         obj.nDoF = 3*obj.totNodMaster+3*obj.totNodSlave+3*obj.nMult;
         obj.E = obj.meshGlue.interfaces.mortar.getMortarOperator(3); % mortar operator
      end


      function computeContactMatrices(obj)
         % get and normalize nodal normals
         n_elems = obj.meshGlue.interfaces(1).elemNormal;
         area_nod = sqrt(sum(obj.meshGlue.interfaces(1).nodeNormal.^2,2)); % nodal area
         n_nodes = obj.meshGlue.interfaces(1).nodeNormal./area_nod;
         % compute mortar matrices within one simulation loop use radial
         % basis functions to evaluate matrix M
         mortar = obj.meshGlue.interfaces(1).mortar; % only one interface defined
         %
         nGP = obj.meshGlue.interfaces(1).nG; % numb. of gauss points for element based integration
         c_ns = 0;  % counter for GP not projected
         % set Gauss class
         gS = Gauss(mortar.slaveCellType,nGP,2); % gauss class for slave integration
         elemSlave = getElem(mortar,gS,'slave');
         [imVec,jmVec,Mgvec,Mtvec,Mnvec] = deal(zeros(nnz(mortar.elemConnectivity)*mortar.nNmaster^2,1));
         [idVec,jsVec,Dgvec,Dtvec,Dnvec] = deal(zeros(nnz(mortar.elemConnectivity)*mortar.nNmaster^2,1));
         [iMapVec,jMapVec,aMapVec] = deal(zeros(nnz(mortar.elemConnectivity),1));

         % Perform interpolation on the master side (computing weights
         % and interpolation coordinates)
         % Loop trough slave elements
         cs = 0; % slave matrix entry counter
         cm = 0; % master matrix entry counter
         cl = 0; % multiplier matrix entry counter
         cmap = 0; % shared are map matrix entry counter
         nEl = mortar.nElSlave;
         for is = 1:nEl
            %Compute Slave quantities
            dJWeighed = elemSlave.getDerBasisFAndDet(is,3); % Weighted Jacobian
            %get Gauss Points position in the real space
            ptsGauss = getGPointsLocation(elemSlave,is);
            nSlave = mortar.intSlave.surfaces(is,:);
            n_el = n_elems(is,:);
            n_nod = (n_nodes(nSlave,:))';
            n_nod = n_nod(:);
            Rloc = getRotationMatrix(obj,n_el);
            %Rloc = eye(3*obj.nNslave);
            %A = diag(repelem(area_nod(nSlave),3));
            NSlave = getBasisFinGPoints(elemSlave); % Get slave basis functions
            master_elems = find(mortar.elemConnectivity(:,is));
            for im = master_elems'
               nMaster = mortar.intMaster.surfaces(im,:);
               [NMaster,id] = mortar.getMasterBasis(im,ptsGauss); % compute interpolated master basis function \Pi(Nm)
               if any(id)
                  % element-based integration
                  % prepare provisional 3D matrices
                  Nm = obj.dispSP(NMaster(id,:));
                  Ns = obj.dispSP(NSlave(id,:));
                  Nmult = ones(size(NSlave(id,:),1),1);
                  Nmult = obj.dispSP(Nmult);
                  % global mortar matrices (global frame)
                  Dgtmp = pagemtimes(Nmult,'transpose',Ns,'none');
                  Dgtmp = Dgtmp.*reshape(dJWeighed(id),1,1,[]);
                  Dgloc = Rloc'*sum(Dgtmp,3);
                  Mgtmp = pagemtimes(Nmult,'transpose',Nm,'none');
                  Mgtmp = Mgtmp.*reshape(dJWeighed(id),1,1,[]);
                  Mgloc = Rloc'*sum(Mgtmp,3);
                  % normal mortar matrices (global frame)
                  Nn = pagemtimes(Ns,n_nod);
                  Dntmp = pagemtimes(Nmult(1,:,:),'transpose',pagemtimes(Nn,'transpose',Ns,'none'),'none');
                  Dntmp = Dntmp.*reshape(dJWeighed(id),1,1,[]);
                  Dnloc = -Rloc'*sum(Dntmp,3);
                  Mntmp = pagemtimes(Nmult(1,:,:),'transpose',pagemtimes(Nn,'transpose',Nm,'none'),'none');
                  Mntmp = Mntmp.*reshape(dJWeighed(id),1,1,[]);
                  Mnloc = -Rloc'*sum(Mntmp,3);
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
                  dof_mult = get_dof(is);
                  [jjM,iiM] = meshgrid(dof_master,dof_mult);
                  [jjD,iiD] = meshgrid(dof_slave,dof_mult);
                  [jjL,iiL] = meshgrid(dof_mult,dof_mult);
                  nm = numel(Mgloc);
                  ns = numel(Dgloc);
                  nl = numel(Lloc);
                  imVec(cm+1:cm+nm) = iiM(:); jmVec(cm+1:cm+nm) = jjM(:);
                  idVec(cs+1:cs+ns) = iiD(:); jsVec(cs+1:cs+ns) = jjD(:);
                  ilVec(cl+1:cl+nl) = iiL(:); jlVec(cl+1:cl+nl) = jjL(:);
                  Mgvec(cm+1:cm+nm) = Mgloc(:);
                  Dgvec(cs+1:cs+ns) = Dgloc(:);
                  Mnvec(cm+1:cm+nm) = Mnloc(:);
                  Dnvec(cs+1:cs+ns) = Dnloc(:);
                  Mtvec(cm+1:cm+nm) = Mtloc(:);
                  Dtvec(cs+1:cs+ns) = Dtloc(:);
                  Lvec(cl+1:cl+nl)  = Lloc(:);
                  % sort out Points already projected
                  iMapVec(cmap+1) = im;
                  jMapVec(cmap+1) = is;
                  aMapVec(cmap+1) = sum(dJWeighed(id));
                  dJWeighed = dJWeighed(~id);
                  ptsGauss = ptsGauss(~id,:);
                  NSlave = NSlave(~id,:);
                  cs = cs+ns;
                  cm = cm+nm;
                  cl = cl+nl;
                  cmap = cmap+1;
               end
            end
            if ~all(id)
               fprintf('GP not sorted for slave elem %i \n',is);
               c_ns = c_ns + 1;
            end
         end
         imVec = imVec(1:cm); jmVec = jmVec(1:cm);
         idVec = idVec(1:cs); jsVec = jsVec(1:cs);
         Mgvec = Mgvec(1:cm); Mnvec = Mnvec(1:cm); Mtvec = Mtvec(1:cm);
         Dgvec = Dgvec(1:cm); Dnvec = Dnvec(1:cm); Dtvec = Dtvec(1:cm);
         iMapVec = iMapVec(1:cmap); jMapVec = jMapVec(1:cmap);
         aMapVec = aMapVec(1:cmap);
         Lvec = Lvec(1:cl);
         obj.Mg = sparse(imVec,jmVec,Mgvec,3*nEl,3*obj.nM);
         obj.Mn = sparse(imVec,jmVec,Mnvec,3*nEl,3*obj.nM);
         obj.Mt = sparse(imVec,jmVec,Mtvec,3*nEl,3*obj.nM);
         obj.Dg = sparse(idVec,jsVec,Dgvec,3*nEl,3*obj.nS);
         obj.Dn = sparse(idVec,jsVec,Dnvec,3*nEl,3*obj.nS);
         obj.Dt = sparse(idVec,jsVec,Dtvec,3*nEl,3*obj.nS);
         obj.L = sparse(ilVec,jlVec,Lvec,3*nEl,3*nEl);
         obj.areaMap = sparse(iMapVec,jMapVec,aMapVec,...
            mortar.nElMaster,mortar.nElSlave);
         obj.areaMap = obj.areaMap./sum(obj.areaMap,2);
         computeStabilizationMatrix(obj);
         % eliminate small numerical quantities in obj.Dg
         %obj.Dg(abs(obj.Dg)<1e-12) = 0;
      end

      function R = getRotationMatrix(obj,nloc)
         % n_el: local nodal normal array
         nNode = round(numel(nloc)/3); 
         R = zeros(3*nNode,3*nNode);
         for i = 0:nNode-1
            n = nloc(3*i+1:3*i+3);
            [Rn,~] = qr(n');
            v = Rn'*n';
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

      function R = getElemRotationMatrix(obj,k)
         % n_el: local nodal normal array
         n = obj.meshGlue.interfaces.elemNormal(k,:);
         n = (n/norm(n,2))';
         [R,~] = qr(n);
         v = R'*n;
         if v(1)>0
            R(:,1) = -R(:,1);
         end
      end
      

      function Rloc = getGlobalRotationMatrix(obj)
         n_el = obj.meshGlue.interfaces(1).elemNormal;
         % n_el: local nodal normal array
         nEntry = 9;
         Rvec = zeros(nEntry*obj.nMult,1);
         iivec = zeros(nEntry*obj.nMult,1);
         jjvec = zeros(nEntry*obj.nMult,1);
         l1 = 0;
         for i = 0:obj.nMult-1
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
         Rloc = sparse(iivec,jjvec,Rvec,3*obj.nMult,3*obj.nMult);
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
         %Lmat = obj.L./sum(obj.L,2);
         %mult = Lmat*mult; % compute variationally consistent traction
         for i = activeSet.curr.slip'
            dofSlip = get_dof(i);
            gapNorm = norm(obj.g_T(dofSlip),2);
            if gapNorm > obj.tolGap
               % principle of maximum plastic dissipation
               tL = repelem(tau_max(obj,i,mult),3,1).*obj.g_T(dofSlip)/gapNorm;
               tL = obj.getNodeRotationMatrix(i)'*tL; % traction in local coords
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
         % compute variationally consistent normal stress
%          Lmat = obj.L./sum(obj.L,2);
%          multipliers = Lmat*multipliers;
         s_n = multipliers(3*nList-2);
         tau = obj.coes - tan(deg2rad(obj.phi))*s_n;
      end

      function computeNodalGap(obj,state,dofMap,mult,varargin)
         % compute gap in global coordinates
         % compute time difference of the tangential component of nodal gap
         % \Delta_n g_T (in global coordinates)
         % Use mortar operator E to map master nodes to slave side
         usCurr = state(obj.tagSlave).dispCurr(dofMap.nodSlave);
         usConv = state(obj.tagSlave).dispConv(dofMap.nodSlave);
         umCurr = state(obj.tagMaster).dispCurr(dofMap.nodMaster);
         umConv = state(obj.tagMaster).dispConv(dofMap.nodMaster);
         % compute weighted nodal gap
         % stabilization contribution to the gap
         stabGap = obj.stabMat*mult;
         if ~isempty(varargin)
            dofStick = varargin{1};
            l = 1:size(obj.stabMat,1);
            id = ~ismember(l,dofStick);
            % careful: we are removing stabilization contribution to open
            % nodes
            stabGap(id) = 0;
         end
         obj.gap = (obj.Dg*usCurr - obj.Mg*umCurr - stabGap)./sum(obj.Dg,2);
         slip = (obj.Dg*(usCurr-usConv) - obj.Mg*(umCurr-umConv) - stabGap)./sum(obj.Dg,2);
         if isempty(obj.gapOld)
            obj.gapOld = obj.gap;
         end
         obj.g_T = slip; % global nodal gap
         n = obj.meshGlue.interfaces.elemNormal;
         % compute tangential component
         obj.g_N = zeros(obj.nMult,1);
         l1 = 0;
         for i = 1:obj.nMult
            % get node normal
            n_i = n(i,:);
            T = eye(3) - n_i'*n_i; % tangential projection matrix
            obj.g_T(l1+1:l1+3) = T*obj.g_T(l1+1:l1+3);
            obj.g_N(i) = -n_i*obj.gap(l1+1:l1+3);
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

      function [TD,TM,N] = computeConsistencyMatrices(obj,activeSet,mult,state,stateOld,dofMap)
         % compute the non linear consistency matrices by directly evaluating
         % the gap at the gauss points of the slave side
         usDiff = state(obj.tagSlave).dispCurr(dofMap.nodSlave)-stateOld(obj.tagSlave).dispCurr(dofMap.nodSlave);
         umDiff = state(obj.tagMaster).dispCurr(dofMap.nodMaster)-stateOld(obj.tagMaster).dispCurr(dofMap.nodMaster);
         n_elems = obj.meshGlue.interfaces(1).elemNormal;
         area_nod = sqrt(sum(obj.meshGlue.interfaces(1).nodeNormal.^2,2)); % nodal area
         n_nodes = obj.meshGlue.interfaces(1).nodeNormal./area_nod;
         mortar = obj.meshGlue.interfaces(1).mortar; % only one interface defined
         %nGP = obj.meshGlue.interfaces(1).nG; % numb. of gauss points for element based integration (use this for all matrices)
         %Mdetect = zeros(mortar.nElMaster,mortar.nElSlave);
         % set Gauss class
         gS = Gauss(mortar.slaveCellType,4,2); % gauss class for slave integration
         elemSlave = getElem(mortar,gS,'slave');
         [imVec,jmVec,Mtvec] = deal(zeros(nnz(mortar.elemConnectivity)*mortar.nNmaster^2,1));
         [isVec,jsVec,Dtvec] = deal(zeros(nnz(mortar.elemConnectivity)*mortar.nNslave^2,1));
         [inVec,jnVec,Nvec] = deal(zeros(mortar.nElSlave*mortar.nNslave^2,1));
         cn = 0; % counter for normal matrix
         cm = 0;
         cs = 0;
         for is = activeSet.curr.slip'
            %Compute Slave quantities
            if norm(obj.g_T(get_dof(is))) < obj.tolGap
               % gap too small to compute consistency matrices
               continue
            end
            dJWeighed = elemSlave.getDerBasisFAndDet(is,3); % Weighted Jacobian
            %get Gauss Points position in the real space
            ptsGauss = getGPointsLocation(elemSlave,is);
            nSlave = mortar.intSlave.surfaces(is,:);
            n_el = n_elems(is,:);
            n_nod = (n_nodes(nSlave,:))';
            n_nod = n_nod(:);
            Rloc = getRotationMatrix(obj,n_el);
            NSlave = getBasisFinGPoints(elemSlave); % Get slave basis functions
            NSlaveMult = ones(size(NSlave,1),1);
            % compute slave only matrices
            Ns = obj.dispSP(NSlave);
            Nmult = obj.dispSP(NSlaveMult);
            gapGP = computeGapInGP(obj,mortar,Ns,usDiff,umDiff,is,ptsGauss);
            Nn = pagemtimes(Ns,n_nod);
            % tangential mortar matrices (global frame)
            Nt = eye(3) - pagemtimes(Nn,'none',Nn,'transpose');
            gtGP = pagemtimes(Nt,gapGP);
            dtdgt = computeDerTracGap(obj,is,mult,gtGP); % 3D matrix with NL derivatives of tractions
            dtdtn = computeDerTracTn(obj,gtGP); % 3D matrix with NL derivatives of tractions
            Ntmp = pagemtimes(Nmult,'transpose',pagemtimes(dtdtn,Nmult(1,:,:)),'none'); % in global coords
            Ntmp =  Ntmp.*reshape(dJWeighed,1,1,[]);
            Nloc = Rloc'*sum(Ntmp,3);
            dof_slave = get_dof(nSlave);
            dof_mult = get_dof(is);
            [jjS,iiS] = meshgrid(dof_slave,dof_mult);
            [jjMult,iiMult] = meshgrid(dof_mult,dof_mult);
            nN = numel(Nloc);
            Nvec(cn+1:cn+nN) = Nloc(:);
            jnVec(cn+1:cn+nN) = jjMult(:); inVec(cn+1:cn+nN) = iiMult(:);
            cn = cn + nN;
            % compute cross-grid matrices
            master_elems = find(mortar.elemConnectivity(:,is)); 
            for im = master_elems'
               nMaster = mortar.intMaster.surfaces(im,:);
               [NMaster,id] = mortar.getMasterBasis(im,ptsGauss); % compute interpolated master basis function \Pi(Nm)
               if any(id)
                  % element-based integration
                  % prepare provisional 3D matrices
                  Nm = obj.dispSP(NMaster(id,:));
                  Ns = obj.dispSP(NSlave(id,:));
                  Nmult = obj.dispSP(NSlaveMult(id,:));
                  dtdgtTmp = dtdgt(:,:,id);
                  % normal mortar matrices (global frame)
                  Nn = pagemtimes(Ns,n_nod);
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
               fprintf('GP not sorted for slave elem %i \n',is);
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


      function rhs = computeRhsLimitTraction2(obj,mult,state,stateOld,activeSet,dofMap)
         % compute rhs related to limit tracion: <mu,t*>_slip
         % gap is comptuted on each gauss point using exact Pi operator
         usDiff = state(obj.tagSlave).dispCurr(dofMap.nodSlave)-stateOld(obj.tagSlave).dispCurr(dofMap.nodSlave);
         umDiff = state(obj.tagMaster).dispCurr(dofMap.nodMaster)-stateOld(obj.tagMaster).dispCurr(dofMap.nodMaster);
         mortar = obj.meshGlue.interfaces(1).mortar; % only one interface defined
         %nGP = obj.meshGlue.interfaces(1).nG; % numb. of gauss points for element based integration (use this for all matrices)
         %Mdetect = zeros(mortar.nElMaster,mortar.nElSlave);
         area_nod = sqrt(sum(obj.meshGlue.interfaces(1).nodeNormal.^2,2)); % nodal area
         n_nodes = obj.meshGlue.interfaces(1).nodeNormal./area_nod; % normalized nodal normals
         n_elems = obj.meshGlue.interfaces(1).elemNormal;
         % set Gauss class
         gS = Gauss(mortar.slaveCellType,4,2); % gauss class for slave integration
         elemSlave = getElem(mortar,gS,'slave');
         rhs = zeros(length(mult),1);
         tLim = obj.computeLimitTraction(activeSet,mult);
         k = 0;
         for j = activeSet.curr.slip'
            % Compute Slave quantities
            % if norm(obj.g_T(get_dof(j))) < obj.tolGap
            %    % gap too small to compute consistency matrices
            %    continue
            % end
            nSlave = mortar.intSlave.surfaces(j,:);
            n_el = n_elems(j,:);
            dofSlave = get_dof(j);
            t = tLim(get_dof(k+1));
            %sn = mult(dofSlave(1:3:end));
            tauLimit = tau_max(obj,j,mult);
            dJWeighed = elemSlave.getDerBasisFAndDet(j,3); % Weighted Jacobian
            %get Gauss Points position in the real space
            ptsGauss = getGPointsLocation(elemSlave,j);
            NSlave = getBasisFinGPoints(elemSlave); % Get slave basis functions
            NSlaveMult = ones(size(NSlave,1),1);
            Rloc = getRotationMatrix(obj,n_el);
            % prepare 3D matrices for Gauss integration    
            if norm(obj.g_T(get_dof(j))) < obj.tolGap
               Nmult = obj.dispSP(NSlaveMult);
               tSlip = pagemtimes(Nmult,t);
            else
               Ns = obj.dispSP(NSlave);
               Nmult = obj.dispSP(NSlaveMult);
               n_nod = (n_nodes(nSlave,:))';
               n_nod = n_nod(:); % local nodal normal vector
               tauLim = pagemtimes(Nmult(1,1:3:end,:),tauLimit);
               gapGP = computeGapInGP(obj,mortar,Ns,usDiff,umDiff,j,ptsGauss);
               % project gap onto tangential direction
               % normal mortar matrices (global frame)
               Nn = pagemtimes(Ns,n_nod);
               % tangential mortar matrices (global frame)
               Nt = eye(3) - pagemtimes(Nn,'none',Nn,'transpose');
               gtGP = pagemtimes(Nt,gapGP);
               norm_gt = pagenorm(gtGP);
               tSlip = pagemtimes(tauLim,pagemtimes(gtGP,1/norm_gt)); 
            end
            rhsTmp = pagemtimes(Nmult(2:3,:,:),'transpose',tSlip(2:3,:,:),'none');
            rhsTmp = rhsTmp.*reshape(dJWeighed,1,1,[]);
            rhsLoc = sum(rhsTmp,3);
            rhs(dofSlave) = rhs(dofSlave)+Rloc'*rhsLoc;
            k = k+1;
         end
      end

      function computeStabilizationMatrix(obj)
         % pressure - jump stabilization matrix in 3D
         % leverage edge topology of mortar 3D matrix
         mortar = obj.meshGlue.interfaces.mortar;
         elem = Elements(mortar.intSlave,Gauss(12,3,2));
         % list of internal edges
         inner_edges = find(all(mortar.edge2cells,2)); 
         iVec = zeros(6*numel(inner_edges),1);
         jVec = iVec; 
         hVec = iVec;
         c = 0;
         for i = inner_edges'
            es = mortar.edge2cells(i,:);
            n = mortar.edge2nodes(i,:);
            Kslave = zeros(3);
            Kmaster = zeros(3);
            id = [ismember(mortar.intSlave.surfaces(es(1),:),n(1));
               ismember(mortar.intSlave.surfaces(es(1),:),n(2))];
            A1 = elem.quad.findNodeArea(es(1));
            A11 = A1(id(1,:)); A12 = A1(id(2,:));
            id = [ismember(mortar.intSlave.surfaces(es(2),:),n(1));
               ismember(mortar.intSlave.surfaces(es(2),:),n(2))];
            A2 = elem.quad.findNodeArea(es(2));
            A21 = A2(id(1,:)); A22 = A2(id(2,:));
            A = A11*A21+A12*A22;
            % get stiffness approximation from master and slave side
            for e = es
               Kslave = Kslave + getStiffSlave(obj,e);
            end
            em = unique([find(obj.areaMap(:,es(1))); find(obj.areaMap(:,es(2)))]);
            for e = em'
               Kmaster = Kmaster + getStiffMaster(obj,e,es);
            end
            S = A*(inv(Kslave) + inv(Kmaster));
            Hloc = 10*[S -S; -S S];
            dof = get_dof(es);
            [jjLoc,iiLoc] = meshgrid(dof,dof);
            cc = numel(Hloc);
            % get node area
            iVec(c+1:c+cc) = iiLoc(:);
            jVec(c+1:c+cc) = jjLoc(:);
            hVec(c+1:c+cc) = Hloc(:);
            c = c+cc;
         end
         obj.stabMat = sparse(iVec,jVec,hVec,3*mortar.nElSlave,3*mortar.nElSlave);
      end


      function computeStabilizationMatrixNew(obj)
         % pressure - jump stabilization matrix in 3D
         % leverage edge topology of mortar 3D matrix
         mortar = obj.meshGlue.interfaces.mortar;
         elem = Elements(mortar.intSlave,Gauss(12,3,2));
         % list of internal edges
         inner_edges = find(all(mortar.edge2cells,2));
         iVec = zeros(6*numel(inner_edges),1);
         jVec = iVec;
         hVec = iVec;
         c = 0;
         for i = inner_edges'
            es = mortar.edge2cells(i,:);
            n = mortar.edge2nodes(i,:);
            Kslave = zeros(3);
            Kmaster = zeros(3);
            id = [ismember(mortar.intSlave.surfaces(es(1),:),n(1));
               ismember(mortar.intSlave.surfaces(es(1),:),n(2))];
            A1 = elem.quad.findNodeArea(es(1));
            A11 = A1(id(1,:)); A12 = A1(id(2,:));
            id = [ismember(mortar.intSlave.surfaces(es(2),:),n(1));
               ismember(mortar.intSlave.surfaces(es(2),:),n(2))];
            A2 = elem.quad.findNodeArea(es(2));
            A21 = A2(id(1,:)); A22 = A2(id(2,:));
            A = A11*A21+A12*A22;
            % get stiffness approximation from master and slave side
            for e = es
               Kslave = Kslave + getStiffSlave(obj,e);
            end
            em = unique([find(obj.areaMap(:,es(1))); find(obj.areaMap(:,es(2)))]);
            for e = em'
               Kmaster = Kmaster + getStiffMaster(obj,e,es);
            end
            S = A*(inv(Kslave) + inv(Kmaster));
            Hloc = [S -S; -S S];
            dof = get_dof(es);
            [jjLoc,iiLoc] = meshgrid(dof,dof);
            cc = numel(Hloc);
            % get node area
            iVec(c+1:c+cc) = iiLoc(:);
            jVec(c+1:c+cc) = jjLoc(:);
            hVec(c+1:c+cc) = Hloc(:);
            c = c+cc;
         end
         obj.stabMat = sparse(iVec,jVec,hVec,3*mortar.nElSlave,3*mortar.nElSlave);
      end


   end

   methods (Access=private)
      function K = getStiffSlave(obj,el)
         % get approximation of stiffess contribution from slave side to
         % the pressure-jump stabilization parameter
         slaveTag = obj.meshGlue.interfaces.Slave;
         mat = obj.meshGlue.model(slaveTag).Material;
         mortar = obj.meshGlue.interfaces.mortar;
         c = mortar.f2cSlave(el);
         cTag = mortar.mshSlave.cellTag(c);
         youngMod = mat.getMaterial(cTag).ConstLaw.E; % young modulus
         vol = mortar.mshSlave.cellVolume(c);
         % size of the bounding box
         x = mortar.mshSlave.coordinates(mortar.mshSlave.cells(c,:),1);
         y = mortar.mshSlave.coordinates(mortar.mshSlave.cells(c,:),2);
         z = mortar.mshSlave.coordinates(mortar.mshSlave.cells(c,:),3);
         lx = abs(max(x) - min(x));
         ly = abs(max(y) - min(y));
         lz = abs(max(z) - min(z));
         K = diag((youngMod*vol)./([lx;ly;lz]));
      end

      function K = getStiffMaster(obj,el,elSlave)
         % get approximation of stiffess contribution from slave side to
         % the pressure-jump stabilization parameter
         r = sum(obj.areaMap(el,elSlave),"all"); % effective area of master element
         masterTag = obj.meshGlue.interfaces.Master;
         mat = obj.meshGlue.model(masterTag).Material;
         mortar = obj.meshGlue.interfaces.mortar;
         c = mortar.f2cMaster(el);
         cTag = mortar.mshSlave.cellTag(c);
         youngMod = mat.getMaterial(cTag).ConstLaw.E; % young modulus
         vol = mortar.mshMaster.cellVolume(c);
         % size of the bounding box
         x = mortar.mshMaster.coordinates(mortar.mshMaster.cells(c,:),1);
         y = mortar.mshMaster.coordinates(mortar.mshMaster.cells(c,:),2);
         z = mortar.mshMaster.coordinates(mortar.mshMaster.cells(c,:),3);
         lx = abs(max(x) - min(x));
         ly = abs(max(y) - min(y));
         lz = abs(max(z) - min(z));
         K = diag((r*youngMod*vol)./([lx;ly;lz]));
      end

      function dtdgt = computeDerTracGap(obj,elemId,mult,gapGPs)
         % N: 3D slave side matrix
         % result 3D matrix of size (3*nN)x(3*nN)xnG
         % nodeId: node list of input element
         nG = size(gapGPs,3);
         % gt = obj.g_T(get_dof(nodeId));
         sigma_n = mult(3*elemId-2);
         dtdgt = repmat(zeros(3,3),1,1,nG);
         for i = 1:nG
            % get limit tau at gauss point
            g = gapGPs(:,:,i);
            tauLim = obj.coes - tan(deg2rad(obj.phi))*sigma_n;
            dtdgt(:,:,i) = tauLim*((eye(3)*norm(g)^2 - g*g')/(norm(g))^3);
         end
      end

      function dtdtn = computeDerTracTn(obj,gap)
         nG = size(gap,3);
         % Nslave, Nmult: 3D slave side matrix
         % result 3D matrix of size (3*nN)x(3*nN)xnG
         dtdtn = repmat(zeros(3,1),1,1,nG);
         for i = 1:nG
            g = gap(:,:,i);
            dtdtn(:,:,i) = -tan(deg2rad(obj.phi))*(g/norm(g));
         end
      end

      function gapGP = computeGapInGP(obj,mortar,Ns,usDiff,umDiff,j,ptsGauss)
          % us and um are the displacement increment w.r.t previously
          % converged time step
         nSlave = mortar.intSlave.surfaces(j,:); % slave node id
         us = pagemtimes(Ns,usDiff(get_dof(nSlave)));
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
               um(:,:,idM(id)) = pagemtimes(Nm,umDiff(get_dof(nMaster)));
            end
            ptsGauss = ptsGauss(~id,:);
            % index of non projected gp from previous projection
            idM = idM(~id);
         end
         gapGP = us - um;
      end
   end
end
