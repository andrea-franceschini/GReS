classdef mortarFaults_old < handle
   properties
      simParameters
      t = 0
      nDom = 2
      nInterf = 1
      tStep = 0
      dt
      state
      stateOld
      models    
      maxActiveSetIters = 4 % maximum number of active set iterations
      meshGlue    % instance of the meshGlue class (contains mortar algorithms)
      multipliers % store current and previous contact tractions
      dofMap      % group DoF based on master/slave/state
      activeSet   % group nodes based on active set state
      gap         % g
      g_T         % \Delta_n g_T
      g_N         % n \cdot g
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
      totNodMaster % total number of master domain nodes
      totNodSlave  % total number of slave domain nodes
      nDoF
      E           % mortar operator
      R           % global rotation matrix
      tolGap = 1e-7; % Gap tolerance to avoid instabilities [m]
      iniMult; % initial multiplier vector
   end

   methods
      function obj = mortarFaults_old(simParam,meshGlue,coes,phi)
         obj.coes = coes;
         obj.phi = phi;
         obj.setParams(simParam,meshGlue);
         obj.simulationLoop();
      end

      function simulationLoop(obj)
         %
         obj.dt = obj.simParameters.dtIni;
         obj.dt; 

         % compute mechanics matrices (linear)
         for i = 1:obj.nDom
            getPoro(obj.models(i).Discretizer).computeMat(obj.state(i),obj.dt);
         end

         flConv = true; %convergence flag

         % Compute contact matrices (only at initialization step)
         computeContactMatrices(obj);

         % Initialize nodal gap
         computeNodalGap(obj);
         
         % TIME LOOP
         while obj.t < obj.simParameters.tMax
            % Update the simulation time and time step ID
            absTol = obj.simParameters.absTol;
            obj.tStep = obj.tStep + 1;
            obj.t = obj.t + obj.dt;
            % Apply Dirichlet value to individual domain solutions
            for i = 1:obj.nDom
               % Apply the Dirichlet condition value to the solution vector
               if ~isempty(obj.meshGlue.model(i).BoundaryConditions)
                  applyDirVal(obj.meshGlue.model(i).ModelType,obj.meshGlue.model(i).BoundaryConditions,...
                     obj.t, obj.state(i));
               end
               %
            end

            if obj.simParameters.verbosity > 0
               fprintf('\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,obj.dt);
               fprintf('-----------------------------------------------------------\n');
            end

            itAS = 0;
            flagActiveSet = true;

            % Active set loop
            while (flagActiveSet) && (itAS < obj.maxActiveSetIters)
               fprintf('Active set iteration n. %i \n',itAS)

               obj.updateDoFMap;

               % Compute rhs terms (no BCS)
               rhs = computeRhs(obj);

               % assemble global jacobian
               J = computeJacobian(obj);

               % apply BCs to global system
               [J,rhs] = applyBC_faults(obj,J,rhs,obj.t);


               % compute base rhs norm for relative convergence
               % Norm of unbalanced forces
               rhsNorm = norm(rhs(1:3*(obj.totNodMaster+obj.totNodSlave)),2);
               itNR = 0;

               % NR printing
               if obj.simParameters.verbosity > 1
                   fprintf('Iter     ||rhs||\n');
               end
               % Relative residual
               tolWeigh = obj.simParameters.relTol*rhsNorm;


               if obj.simParameters.verbosity > 1
                   fprintf('0     %e\n',rhsNorm);
               end

               % Newton Raphson loop
               while ((rhsNorm > tolWeigh) && (itNR < obj.simParameters.itMaxNR) ...
                       && (rhsNorm > absTol)) || itNR == 0
                   itNR = itNR + 1;
                   %
                   du = J\(-rhs);

                   % update solution fields (disp, multipliers)
                   updateStateFaults(obj,du);

                   % printe state (even if not converged)
                   for i = 1:obj.nDom
                       %obj.state(i).advanceState();
                       printState(obj.models(i).OutState,obj.state(i));
                       obj.models(i).OutState.finalize();
                       plotFunction(obj.meshGlue.interfaces.mortar.intSlave,'Fault',(reshape(obj.multipliers,3,[]))',...
                           ["sigma_n","tau_1","tau_2"]);
                   end

                   % update nodal gaps
                   computeNodalGap(obj);

                   % recompute rhs
                   rhs = computeRhs(obj);

                   % assemble jacobian
                   J = computeJacobian(obj);

                   % Apply BCs to global system
                   [J,rhs] = applyBC_faults(obj,J,rhs,obj.t);

                   % compute Rhs norm
                   rhsNorm = norm(rhs,2);
                   if obj.simParameters.verbosity > 1
                       fprintf('%d     %e\n',itNR,rhsNorm);
                   end
               end
               %
               % Check NR convergence
               flConv = (rhsNorm < tolWeigh || rhsNorm < absTol);
               if flConv % Convergence
                   for i = 1:obj.nDom
                       obj.state(i).curr.t = obj.t;
                       % Print the solution, if needed
                       obj.state(i).advanceState();
                   end
               if obj.t > obj.simParameters.tMax   % For Steady State
                   for i = 1:obj.nDom
                       printState(obj.models(i).OutState,obj.state(i));
                   end
               else
                   for i = 1:obj.nDom
                       printState(obj.models(i).OutState,obj.state(i).prev,obj.state(i).curr);
                   end
               end
               % Update active set
               flagActiveSet = updateActiveSet(obj);
               end
            end
            %
            % Manage next time step
            delta_t = manageNextTimeStep(obj,delta_t,flConv);
         end
         %
      end

      function setParams(obj,simParam,mG)
         obj.simParameters = simParam;
         obj.models = mG.model;
         obj.meshGlue = mG;
         obj.nNmaster = mG.interfaces.mortar.nNmaster;
         obj.nNslave = mG.interfaces.mortar.nNslave;
         assert(numel(obj.models)==2,'Too many input domains');
         % pattern of basis functions transform
         obj.indN = repmat([1 5 9],1,obj.nNslave)+9*repelem(0:obj.nNslave-1,1,3);
         % save model metrics
         obj.tagMaster = obj.meshGlue.interfaces(1).Master;
         obj.tagSlave = obj.meshGlue.interfaces(1).Slave;
         obj.totNodMaster = obj.models(obj.tagMaster).Grid.topology.nNodes;
         obj.totNodSlave = obj.models(obj.tagSlave).Grid.topology.nNodes;
         obj.nS = length(mG.interfaces.mortar.nodesSlave);
         obj.nM = length(mG.interfaces.mortar.nodesMaster);
         idMaster = obj.meshGlue.interfaces(1).masterSet;
         idSlave = obj.meshGlue.interfaces(1).slaveSet;
         % add lagrange multiplier field to state structure
         obj.multipliers = zeros(3*obj.nS,1);
         obj.gap = zeros(3*obj.nS,1);  % NO INITIAL GAP
         % state struct easier to access (we only work with a curr. state)
         obj.state = [obj.models(1).State;obj.models(2).State];
         obj.stateOld = [copy(obj.models(1).State);copy(obj.models(2).State)];
         % compute global rotation matrix
         obj.R = getGlobalRotationMatrix(obj);
         % initialize stress (not flexible at all)
         setInitialStress(obj,'x',-1);
         % initialize active set structure
         obj.activeSet = struct('prev',[],'curr',[]);
         obj.activeSet.prev = struct('stick',(1:obj.nS)','slip',[],'open',[]);
         obj.activeSet.curr = obj.activeSet.prev;
         % build dof map with initial active set status
         obj.nDoF = 3*obj.totNodMaster+3*obj.totNodSlave+3*numel(idSlave);
         buildDoFMap(obj,idMaster,idSlave);
         obj.E = obj.meshGlue.interfaces.mortar.getMortarOperator(3); % mortar operator
      end

      function [dt] = manageNextTimeStep(obj,dt,flConv)
         if ~flConv   % Perform backstep
            for i = 1:obj.nDom
               transferState(obj.state(i).curr,obj.state(i).prev);
            end
            obj.t = obj.t - obj.dt;
            obj.tStep = obj.tStep - 1;
            dt = dt/obj.simParameters.divFac;
            obj.dt = obj.dt/obj.simParameters.divFac;  % Time increment chop
            if min(dt,obj.dt) < obj.simParameters.dtMin
               if obj.simParameters.goOnBackstep == 1
                  flConv = 1;
               elseif obj.simParameters.goOnBackstep == 0
                  error('Minimum time step reached')
               end
            elseif obj.simParameters.verbosity > 0
               fprintf('\n %s \n','BACKSTEP');
            end
         end
         if flConv % Go on if converged
            obj.dt = min([obj.dt*obj.simParameters.multFac,obj.simParameters.dtMax]);
            obj.dt = max([obj.dt obj.simParameters.dtMin]);
            for i = 1:obj.nDom
               transferState(obj.state(i).curr,obj.state(i).prev);
            end
            %
            if ((obj.t + obj.dt) > obj.simParameters.tMax)
               obj.dt = obj.simParameters.tMax - obj.t;
            end
         end
      end

      function updateStateFaults(obj,du)
         % update solution vectors and derived quantities (strains)
         obj.state(obj.tagMaster).dispCurr = obj.state(obj.tagMaster).dispCurr + ...
            du(obj.dofMap.master);
         obj.state(obj.tagSlave).dispCurr = obj.state(obj.tagSlave).dispCurr + ...
            du(obj.dofMap.slave);
         obj.multipliers = obj.multipliers + du(obj.dofMap.slave(end)+1:end);
         %
         % Update stress
         l = 0;
         for i = 1:obj.nDom
            du = obj.state(i).dispCurr - obj.state(i).dispConv;
            for el=1:obj.models(i).Grid.topology.nCells
               dof = getDoFID(obj.models(i).Grid.topology,el);
               switch obj.models(i).Grid.topology.cellVTKType(el)
                  case 10 % Tetra
                     N = getDerBasisF(obj.models(i).Grid.cells.tetra,el);
                     B = zeros(6,4*obj.models(i).Grid.topology.nDim);
                     B(obj.models(i).Grid.cells.indB(1:36,2)) = ...
                        N(obj.models(i).Grid.cells.indB(1:36,1));
                     obj.state(i).curr.strain(l+1,:) = (B*du(dof))';
                     l = l + 1;
                  case 12 % Hexa
                     N = getDerBasisFAndDet(obj.models(i).Grid.cells.hexa,el,2);
                     B = zeros(6,8*obj.models(i).Grid.topology.nDim,obj.models(i).Gauss.nNode);
                     B(obj.models(i).Grid.cells.indB(:,2)) = N(obj.models(i).Grid.cells.indB(:,1));
                     obj.state(i).curr.strain((l+1):(l+obj.models(i).Gauss.nNode),:) = ...
                        reshape(pagemtimes(B,du(dof)),6,obj.models(i).Gauss.nNode)';
                     l = l + obj.models(i).Gauss.nNode;
               end
            end
         end
      end



      function [t, dt] = updateTime(obj,conv,dt)
         tMax = obj.simParameters.tMax;
         for i = 1:obj.nDom
            if obj.meshGlue.model(i).OutState.modTime
               tmp = find(obj.t<obj.meshGlue.model(i).outState.timeList(),1,'first');
               if ~conv
                  t = min([obj.t + obj.dt, obj.t + dt, obj.meshGlue.model(i).OutState.timeList(tmp)]);
               else
                  t = min([obj.t + obj.dt, obj.meshGlue.model(i).OutState.timeList(tmp)]);
               end
            else
               t = obj.t + obj.dt;
            end
            if t > tMax
               t = tMax;
            end
         end
         dt = t - obj.t;
      end


      function computeContactMatrices(obj)
         % get and normalize nodal normals
         area_nod = sqrt(sum(obj.meshGlue.interfaces(1).nodeNormal.^2,2)); % nodal area
         n = obj.meshGlue.interfaces(1).nodeNormal./area_nod;
         % compute mortar matrices within one simulation loop use radial
         % basis functions to evaluate matrix M
         mortar = obj.meshGlue.interfaces(1).mortar; % only one interface defined
         nGP = 4; % numb. of gauss points for element based integration
         nInt = 4; % numb. of interpolation points for RBF intrpolation
         tol = 1e-3;
         type = 'gauss';
         mult_type = 'dual';
         c_ns = 0;  % counter for GP not projected
         %Mdetect = zeros(mortar.nElMaster,mortar.nElSlave);
         % set Gauss class
         gM = Gauss(mortar.masterCellType,3,2); % gauss class for Master element interpolation
         gS = Gauss(mortar.slaveCellType,nGP,2); % gauss class for slave integration
         elemMaster = getElem(mortar,gM,'master');
         elemSlave = getElem(mortar,gS,'slave');
         [imVec,jmVec,Mgvec,Mtvec,Mnvec] = deal(zeros(nnz(mortar.elemConnectivity)*mortar.nNmaster^2,1));
         [isVec,jsVec,Dgvec,Dtvec,Dnvec] = deal(zeros(nnz(mortar.elemConnectivity)*mortar.nNmaster^2,1));
         % Perform interpolation on the master side (computing weights
         % and interpolation coordinates)
         [wFMat,w1Mat,ptsIntMat] = mortar.getWeights('master',nInt,elemMaster,type);
         %[wFMatS,w1MatS,ptsIntMatS] = getWeights(obj,'slave',nInt,elemSlave,type);
         % Interpolation for support detection
         %[wFSupp,w1Supp] = getSuppWeight(obj);
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
            Rloc = -getRotationMatrix(obj,n_el);
            %Rloc = eye(3*obj.nNslave);
            %A = diag(repelem(area_nod(nSlave),3));
            NSlave = getBasisFinGPoints(elemSlave); % Slave basis functions
            switch mult_type
               case 'standard'
                  NSlaveMult = NSlave; % Slave basis functions
               case 'dual'
                  NSlaveMult = mortar.computeDualBasisF(NSlave,dJWeighed);
            end
            master_elems = find(mortar.elemConnectivity(:,j));
            for jm = master_elems'
               nMaster = mortar.intMaster.surfaces(jm,:);
               ptsInt = ptsIntMat(:,repNum(3,jm));
               [fiNM,id1] = mortar.computeRBFfiNM(ptsInt,ptsGauss,type);
               switch mortar.degree
                  case 1
                     NMaster = (fiNM*wFMat(:,repNum(mortar.nNmaster,jm)))./(fiNM*w1Mat(:,jm));
                     Nsupp = NMaster(:,[1 2 3]);
                  case 2
                     Ntmp = (fiNM*wFMat(:,repNum(mortar.nNmaster+2,jm)))./(fiNM*w1Mat(:,jm));
                     NMaster = Ntmp(:,1:mortar.nNmaster);
                     Nsupp = Ntmp(:,[end-1 end]);
               end
               % automatically detect supports computing interpolant
               id = all([Nsupp >= 0-tol id1],2);
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
                  Nt = repmat(eye(3),[1,1,[]]) - pagemtimes(Nn,'none',Nn,'transpose');
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
                  dof_master = obj.get_dof(nMaster);
                  dof_slave = obj.get_dof(nSlave);
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
         obj.Dg(abs(obj.Dg)<1e-12) = 0;
      end

      function R = getRotationMatrix(obj,n_el)
         % n_el: local nodal normal array
         R = zeros(3*obj.nNslave,3*obj.nNslave);
         for i = 0:obj.nNslave-1
            n = n_el(3*i+1:3*i+3);
            if all(abs(n) ~= [1;0;0])
               u = [1; 0; 0];
            else
               u = [0; 1; 0];
            end
            m1 = cross(n,u)./norm(cross(n,u),2);
            m2 = cross(n,m1);
            R(3*i+1:3*i+3,3*i+1:3*i+3) = [n m1 m2];
         end
      end

      function R = getNodeRotationMatrix(obj,k)
         % n_el: local nodal normal array
         n = obj.meshGlue.interfaces.nodeNormal(k,:);
         n = (n/norm(n,2))';
         if all(abs(n) ~= [1;0;0])
            u = [1; 0; 0];
         else
            u = [0; 1; 0];
         end
         m1 = cross(n,u)./norm(cross(n,u),2);
         m2 = cross(n,m1);
         R = [n m1 m2];
      end
      

      function R = getGlobalRotationMatrix(obj)
         area_nod = sqrt(sum(obj.meshGlue.interfaces(1).nodeNormal.^2,2));
         n_el = obj.meshGlue.interfaces(1).nodeNormal./area_nod; % normalized nodal area
         % n_el: local nodal normal array
         nEntry = 9;
         Rvec = zeros(nEntry*obj.nS,1);
         iivec = zeros(nEntry*obj.nS,1);
         jjvec = zeros(nEntry*obj.nS,1);
         l1 = 0;
         for i = 0:obj.nS-1
            n = n_el(i+1,:);
            if all(abs(n) ~= [1;0;0])
               u = [1; 0; 0];
            else
               u = [0; 1; 0];
            end
            m1 = cross(n,u)./norm(cross(n,u),2);
            m2 = cross(n,m1);
            [iiLoc,jjLoc] = meshgrid(3*i+1:3*i+3,3*i+1:3*i+3);
            Rloc = [n' m1' m2'];
            iivec(l1+1:l1+nEntry) = iiLoc(:);
            jjvec(l1+1:l1+nEntry) = jjLoc(:);
            Rvec(l1+1:l1+nEntry) = Rloc(:);
            l1 = l1+nEntry;
         end
         R = sparse(iivec,jjvec,Rvec,3*obj.nS,3*obj.nS);
      end

      function Nout = dispSP(obj,Nin)
         % reshape basis functions to obtain displacement shape function
         % input: nG x nN matrix
         % output: 3 x nN x nG
         s2 = 3*obj.nNslave;
         s3 = size(Nin,1);
         N = zeros(3*s2,s3);
         Nin = repelem(Nin',3,1);
         N(obj.indN,:) = Nin;
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
               dofs = obj.get_dof(dofIn);
            case 'slave'
               if isempty(dofIn)
                  dofIn = 1:obj.totNodSlave;
               end
               dofs = 3*obj.totNodMaster+obj.get_dof(dofIn);
            case 'lag'
               if isempty(dofIn)
                  dofs = [];
                  return
               end
               dofs = 3*(obj.totNodMaster+obj.totNodSlave)+obj.get_dof(dofIn);
            otherwise
               error(['Invalide tag string for getContactDofs method \n:' ...
                  'valid inputs are: master, slave, lag']);
         end
      end

      function rhsVal = computeRhsMechanics(obj)
         % Compute rhs of individual domains
         for i = 1:obj.nDom
            % Update tmpState
            getPoro(obj.models(i).Discretizer).computeRhs(obj.state(i));
         end
         % residual standard internal forces
         rhsMaster = getField(obj.models(obj.tagMaster).Discretizer,'Poro').rhs;
         rhsSlave = getField(obj.models(obj.tagSlave).Discretizer,'Poro').rhs;
         rhsVal = [rhsMaster; rhsSlave];
      end
      
      function buildDoFMap(obj,idMaster,idSlave)
         obj.dofMap = struct('master',[],...
            'slave',[],'intMaster',[],'intSlave',[],'stick',[],'slip',[],...
            'open',[],'nodSlave',[],'nodMaster',[]);
         obj.dofMap.master = getContactDofs(obj,'master');
         obj.dofMap.slave = getContactDofs(obj,'slave');
         obj.dofMap.intMaster = getContactDofs(obj,'master',idMaster');
         obj.dofMap.intSlave = getContactDofs(obj,'slave',idSlave');
         obj.dofMap.stick = getContactDofs(obj,'lag',obj.activeSet.curr.stick);
         obj.dofMap.slip = getContactDofs(obj,'lag',obj.activeSet.curr.slip);
         obj.dofMap.open = getContactDofs(obj,'lag',obj.activeSet.curr.open);
         obj.dofMap.nodSlave = obj.get_dof(obj.meshGlue.interfaces.mortar.nodesSlave);
         obj.dofMap.nodMaster = obj.get_dof(obj.meshGlue.interfaces.mortar.nodesMaster); 
      end

      function updateDoFMap(obj)
         obj.dofMap.stick = getContactDofs(obj,'lag',obj.activeSet.curr.stick);
         obj.dofMap.slip = getContactDofs(obj,'lag',obj.activeSet.curr.slip);
         obj.dofMap.open = getContactDofs(obj,'lag',obj.activeSet.curr.open);
      end

      function tLim = computeLimitTraction(obj)
         % TO DO: add check on magnitud gap
         % get limit traction array in local coordinates
         tLim = zeros(3*numel(obj.activeSet.new.slip),1);
         for i = obj.activeSet.new.slip'
            dofSlip = obj.get_dof(i);
            gapNorm = norm(obj.g_T(i),2);
            if gapNorm < obj.tolGap
               tL = zeros(3,1);
            else
               tL = repelem(tau_max(obj,i),3,1).*obj.g_T(dofSlip)/gapNorm; % traction in global coords
               tL = obj.getNodeRotationMatrix(i)'*tL; % traction in local coords
               tL(1) = 0; % set normal component to 0 for rhs computation
            end
            tLim(3*i-2:3*i) = tL;
         end
      end

      function tau = tau_max(obj,i)
         % use the normal component of contact traction 
         tau = obj.coes - tan(deg2rad(obj.phi))*obj.multipliers(3*i-1);
      end

      function computeNodalGap(obj)
         % compute gap in global coordinates 
         % compute time difference of the tangential component of nodal gap
         % \Delta_n g_T (in global coordinates)
         % Use mortar operator E to map master nodes to slave side
         usCurr = obj.state(obj.tagSlave).dispCurr(obj.dofMap.nodSlave);
         umCurr = obj.state(obj.tagMaster).dispCurr(obj.dofMap.nodMaster);
         g_old = obj.gap;
         obj.gap = usCurr - obj.E*umCurr;
         obj.g_T = obj.gap - g_old; % global nodal gap
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

      function rhs = computeRhs(obj)
         % global rhs computation method 
         % rhs contribution of BCs is considered in another method!
         rhs = zeros(obj.nDoF,1);
         % Mechanics internal forces
         rhs([obj.dofMap.master;obj.dofMap.slave]) = rhs([obj.dofMap.master;obj.dofMap.slave])+ computeRhsMechanics(obj);
         % Mesh tying rhs
         rhs([obj.dofMap.intMaster;obj.dofMap.intSlave]) = ...
            rhs([obj.dofMap.intMaster;obj.dofMap.intSlave]) + computeRhsMeshTying(obj);
         % Stick nodes rhs
         rhs(obj.dofMap.stick) = rhs(obj.dofMap.stick) + computeRhsStick(obj);
         % Slip nodes rhs
         rhs(obj.dofMap.slip) = rhs(obj.dofMap.slip) + computeRhsSlip(obj);
         % Open nodes rhs
         rhs(obj.dofMap.open) = rhs(obj.dofMap.open) + computeRhsOpen(obj);
      end

      function J = computeJacobian(obj)
         % dof ordering: master - slave - lagrange
         % interface dofs are not separated by inner dofs numbering
         i = []; j = []; Jvec = []; 
         lagDof = getContactDofs(obj,'lag',1:obj.nS); % complete set of multipliers
         % mechanics block
         Km = getPoro(obj.models(obj.tagMaster).Discretizer).K;
         Ks = getPoro(obj.models(obj.tagSlave).Discretizer).K;
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.master,obj.dofMap.master,Km);
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slave,obj.dofMap.slave,Ks);
         clear Km; clear Ks;
         % mesh tying blocks
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.intMaster,lagDof,-obj.Mg');
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.intSlave,lagDof,obj.Dg');
         % stick blocks
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.stick,obj.dofMap.intMaster,-obj.Mn);
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.stick,obj.dofMap.intSlave,obj.Dn);
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.stick,obj.dofMap.intMaster,-obj.Mt);
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.stick,obj.dofMap.intSlave,obj.Dt);
         % slip blocks
         if any(obj.activeSet.curr.slip)
            slipDofs = obj.get_dof(obj.activeSet.curr.slip);
            % contribution to normal component
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.intMaster,-obj.Mn(slipDofs,:));
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.intSlave,obj.Dn(slipDofs,:));
            % contribution to tangential component
            T = computeDtDgt(obj); % non linear contribution
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.intMaster,-T*obj.Mt(slipDofs,:));
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.intSlave,T*obj.Dt(slipDofs,:));
            N = computeDtDtn(obj); % non linear contribution
            % select only tangential components of mass matrix
            tComp = repmat([false;true;true],numel(obj.activeSet.curr.slip),1);
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.slip,...
               N*obj.L(slipDofs,slipDofs));
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.slip,...
               tComp'.*obj.L(slipDofs,slipDofs).*tComp);
         end
         % open blocks
         if any(obj.activeSet.curr.open)
            openDofs = obj.get_dof(obj.activeSet.curr.open);
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.open,obj.dofMap.open,...
               obj.L(openDofs,openDofs));
         end
         J = sparse(i,j,Jvec,obj.nDoF,obj.nDoF);
      end


      function rhsVal = computeRhsMeshTying(obj)
         % subtract initial stress state
         rhsMaster = -obj.Mg'*(obj.multipliers-obj.iniMult);
         rhsSlave = obj.Dg'*(obj.multipliers-obj.iniMult);
         rhsVal = [rhsMaster;rhsSlave];
      end

      function updateActiveSet(obj)
         % update and compare active set
         % return true if some changes are made
         % return 0 if the active set stays the same
      end
 
      function rhsStick = computeRhsStick(obj)
         usCurr = obj.state(obj.tagSlave).dispCurr(obj.dofMap.nodSlave);
         usConv = obj.state(obj.tagSlave).dispConv(obj.dofMap.nodSlave);
         umCurr = obj.state(obj.tagMaster).dispCurr(obj.dofMap.nodMaster);
         umConv = obj.state(obj.tagMaster).dispConv(obj.dofMap.nodMaster);
         rhsStick = obj.Dn*usCurr - obj.Mn*umCurr + ...
            obj.Dt*(usCurr-usConv) - obj.Mt*(umCurr-umConv);
         rhsStick = rhsStick(obj.get_dof(obj.activeSet.curr.stick));
      end

      function rhsSlip = computeRhsSlip(obj)
         if ~isempty(obj.activeSet.curr.slip)
            dofSlip = get_dof(obj.activeSet.curr.slip);
            tracLim = computeLimitTraction(obj); % limit traction in active nodes
            usCurr = obj.state(obj.tagSlave).dispCurr(obj.dofMap.nodSlave);
            umCurr = obj.state(obj.tagMaster).dispConv(obj.dofMap.nodMaster);
            rhsSlip1 = obj.Dn*usCurr - obj.Mn*umCurr; % normal components
            rhsSlip2 = obj.L(dofSlip,dofSlip)*(obj.multipliers(dofSlip)-tracLim); % tangential components
            rhsSlip = rhsSlip1(dofSlip) + rhsSlip2;
         else
            rhsSlip = [];
         end
      end

      function rhsOpen = computeRhsOpen(obj)
         if ~isempty(obj.activeSet.curr.open)
         rhsOpen = obj.L*obj.multipliers;
         rhsOpen = rhsOpen(obj.get_dof(obj.activeSet.curr.slip));
         else
            rhsOpen = [];
         end
      end

      function Tmat = computeDtDgt(obj)
         % compute derivative of tangential traction w.r.t tangential gap
         % quantities are in global coordinates
         iVec = zeros(9*numel(obj.activeSet.curr.slip),1);
         jVec = zeros(9*numel(obj.activeSet.curr.slip),1);
         Tvec = zeros(9*numel(obj.activeSet.curr.slip),1);
         l = 0;
         for i = obj.activeSet.curr.slip'
            dof = obj.get_dof(i);
            tLim = tau_max(obj,i);
            DgT = obj.g_T(dof);
            normgT = norm(DgT);
            Tloc = tLim*(normgT^2*eye(3)-Dgt*Dgt')/normgT^3; % 3x3 local mat in global coords
            [ii,jj] = meshgrid(dof,dof);
            iVec(l+1:l+9) = ii(:);
            jVec(l+1:l+9) = jj(:);
            Tvec(l+1:l+9) = Tloc(:);
            l = l+9;
         end
         Tmat = sparse(iVec,jVec,Tvec,3*obj.nS,3*obj.nS);
      end

      function Nmat = computeDtDtn(obj)
         % compute derivative of tangential traction w.r.t tangential gap
         % quantities are in local coordinates
         iVec = zeros(2*numel(obj.activeSet.curr.slip),1);
         jVec = zeros(2*numel(obj.activeSet.curr.slip),1);
         Nvec = zeros(2*numel(obj.activeSet.curr.slip),1);
         l = 0;
         for i = obj.activeSet.curr.slip'
            dof = obj.get_dof(i);
            tLim = tau_max(obj,i);
            DgT = obj.getNodeRotationMatrix(i)'*obj.g_T(dof); % local tangential gap
            DgT = DgT([2 3]); % if evrything is ok, the first component is actually 0
            normgT = norm(DgT);
            Nloc = tLim*gT/normgT; % 3x3 local mat in global coords
            iVec(l+1:l+2) = [dof(2);dof(3)];
            jVec(l+1:l+2) = [dof(1);dof(1)];
            Nvec(l+1:l+2) = Nloc(:);
            l = l+2;
         end
         Nmat = sparse(iVec,jVec,Nvec,3*obj.nS,3*obj.nS);
      end

      function setInitialStress(obj,dir,val)
         % dir: define direction of traction enforment (in global reference)
         % 'x' 'y' 'z'
         % val: traction magnitude 
         id = strcmp(["x","y","z"],dir);
         v = val*id;
         % set normal stress for master and slave domains
         id2 = [id false(1,3)];
         obj.state(1).iniStress(:,id2) = val; 
         obj.state(2).iniStress(:,id2) = val; 
         % refer stress to local nodal normal
         v_loc = obj.R*repmat(v',obj.nS,1);
         % make sure contact traction is negative (compression)
         obj.multipliers(:) = -abs(v_loc);
         obj.iniMult = obj.multipliers;
         obj.state(1).conv.stress = obj.state(1).iniStress;
         obj.state(2).conv.stress = obj.state(2).iniStress;
      end

      function [J,rhs] = applyBC_faults(obj,J,rhs,t)
         % Apply Boundary condition to fault mechanics system
         % Block dof indexing is used employing getContactDoF method of mortar
         % faults class
         % Impose BC to the linearized system (Jacobian matrix + RHS)
         for domID = 1:obj.nDom
            model = obj.models(domID).ModelType;
            bound = obj.models(domID).BoundaryConditions;
            keys = bound.db.keys;
            for i = 1 : length(keys)
               dirVal = []; % if stays empty Penalty method is used
               rhsVal = []; % if stays empty Penalty method is used
               cond = bound.getCond(keys{i});
               type = bound.getType(keys{i});
               switch cond
                  case 'NodeBC'
                     bcDofs = bound.getDofs(keys{i});
                     switch bound.getType(keys{i})
                        case 'Neu'
                           rhsVal = - bound.getVals(keys{i}, t);
                     end
                  case 'SurfBC'
                     if isFEMBased(model,'Poro')
                        bcDofs = bound.getDofs(keys{i});
                        switch type
                           case 'Neu'
                              entitiesInfl = bound.getEntitiesInfluence(keys{i});
                              q = bound.getVals(keys{i}, t);
                              rhsVal = - entitiesInfl*q;
                        end
                     end
               end

               % Map single domain dofs to global linear system
               % depending on the domain type
               if domID == obj.tagMaster
                  dof = bcDofs;
               elseif domID == obj.tagSlave
                   bcDofs = bcDofs + 3*obj.totNodMaster;
                    if strcmp(type,'Dir')  % remove constraint on dofs belonging to the interface
                        dof = bcDofs(~ismember(bcDofs,obj.dofMap.intSlave));
                    else
                       dof = bcDofs;
                   end
               end
               % mapping is trivial with only two domains
               switch type
                  case 'Dir' % Dirichlet BC
                     nrows = size(J,1);
                     if isempty(rhsVal) && isempty(dirVal) % FEM Dirichlet BCs
                        vals = zeros(numel(dof),1);
                        [J,rhs] = applyDir(dof,vals,J,rhs);
                     else
                        J(nrows*(dof-1) + dof) = J(nrows*(dof-1) + dof) + dirVal;
                        rhs(dof) = rhs(dof) + rhsVal;
                     end
                  case 'Neu'
                     rhs(dof) = rhs(dof) + rhsVal;
                  otherwise
                     rhs(dof) = rhs(dof) + rhsVal;
               end
            end
         end
      end
   end

   methods (Static)
      function dof = get_dof(nodList)
         % slaveNodes is a row vector
         nodList = reshape(nodList,1,[]);
         dof = repelem(3*nodList,3);
         dof = (dof + repmat(-2:0,1,numel(nodList)))';
      end
   end
end
