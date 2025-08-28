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
         obj.L = sparse(ilVec,jlVec,Lvec,3*nEl,3*obj.nS);
         obj.areaMap = sparse(iMapVec,jMapVec,aMapVec,...
            mortar.nElMaster,mortar.nElSlave);
         obj.areaMap = obj.areaMap./sum(obj.areaMap,2);
         computeStabilizationMatrix(obj);
         % eliminate small numerical quantities in obj.Dg
         %obj.Dg(abs(obj.Dg)<1e-12) = 0;
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
