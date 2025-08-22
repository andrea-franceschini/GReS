classdef NonLinearSolverFaults < handle
   properties
      simParameters
      mortar
      models
      t = 0
      nDom = 2
      nInterf = 1
      tStep = 0
      dt
      state
      stateOld
      maxActiveSetIters = 4 % maximum number of active set iterations
      activeSet
      iniMult; % initial multiplier vector
      currMultipliers % store current and previous contact tractions
      prevMultipliers
      dofMap      % group DoF based on master/slave/state
      dirichNodes
      results
      multType
      convActiveSet
   end

   methods
      function obj = NonLinearSolverFaults(simParam,mG,c,phi,type)
         obj.simParameters = simParam;
         if strcmp(type,'P0')
            obj.mortar = mortarFaultsP0(mG,c,phi);
         else
            obj.mortar = mortarFaults(mG,c,phi);
         end
         obj.multType = type;
         obj.setNonLinearSolverFaults(mG);
         obj.simulationLoop(mG);
      end

      function simulationLoop(obj,mG)
         % sorting nodes and surfaces on symmetry axis
         coordSlave1 = obj.models(2).Grid.topology.coordinates(obj.mortar.idSlave,:);
         coordSlave2 = obj.models(2).Grid.topology.coordinates;
         id1 = find(abs(coordSlave1(:,2)-5)<1e-2);
         id2 = find(all([abs(coordSlave2(:,2)-5)<1e-2 abs(coordSlave2(:,1)-2.5)<1e-2],2));
         [~,idSort1] = sort(coordSlave1(id1,3),'ascend');
         [~,idSort2] = sort(coordSlave2(id2,3),'ascend');
         id1 = id1(idSort1);
         id2 = id2(idSort2);
         if strcmp(obj.multType,'P0')
            id1 = getSymSurf(mG);
         end
           
         % provisional print utilities for faults
         if strcmp(obj.mortar.meshGlue.interfaces.multType,'dual')
            vtkFault = VTKOutput(obj.mortar.meshGlue.interfaces.mortar.intSlave,'Fault_dual');
         elseif strcmp(obj.mortar.meshGlue.interfaces.multType,'standard')
            vtkFault = VTKOutput(obj.mortar.meshGlue.interfaces.mortar.intSlave,'Fault_standard');
         elseif strcmp(obj.mortar.meshGlue.interfaces.multType,'P0')
            vtkFault = VTKOutput(obj.mortar.meshGlue.interfaces.mortar.intSlave,'Fault_P0');
         end
         tListFault = obj.models(1).OutState.timeList;
         printID = 1;
         %
         obj.dt = obj.simParameters.dtIni;

         % compute mechanics matrices (linear)
         for i = 1:obj.nDom
            getPoro(obj.models(i).Discretizer).computeMat(obj.state(i),obj.dt);
         end

         flConv = true; %convergence flag

         % Compute contact matrices (only at initialization step)
         computeContactMatrices(obj.mortar);
         
         % TIME LOOP
         while obj.t < obj.simParameters.tMax
            % handle issue of crack fully sliding

            % Update the simulation time and time step ID
            absTol = obj.simParameters.absTol;
            % Initialize nodal gap
            computeNodalGap(obj.mortar,obj.state,obj.dofMap,obj.currMultipliers,get_dof(obj.activeSet.curr.stick));
            itAS = 0;

            % initialize active set if using nodal multipliers
            obj.initActiveSet();

            %
            flagActiveSet = false;

            obj.mortar.gapOld = obj.mortar.gap; % previous gap vector

            obj.tStep = obj.tStep + 1;
            % if obj.t < 11 && obj.t + obj.dt> 11.01
            %    obj.t = 11;
            %    obj.simParameters.dtMax = 0.25;
            % else
               obj.t = obj.t + obj.dt;
            %end
            if obj.t>13.9
               obj.maxActiveSetIters = 3;
            end
            % update structure for results printing
            initVecConv = zeros(obj.maxActiveSetIters*obj.simParameters.itMaxNR,1);
            initVecTraction = zeros(numel(id1),1);
            initVecDisp = zeros(numel(id2),1);
            obj.results = [obj.results;...
               struct('itNR',initVecConv,'itAS',initVecConv,...
               'rhsNorm',initVecConv,'s_n',initVecTraction,'tauNorm',initVecTraction,'gap',initVecDisp)];

            if (obj.simParameters.verbosity > 0) && (itAS == 0)
               fprintf('\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,obj.dt);
               fprintf('-----------------------------------------------------------\n');
            end


            % Apply the Dirichlet condition value to the solution vector
            applyDirFault(obj);
            k = 0;
            while flagActiveSet || itAS==0
               % Newton Raphson loop
               fprintf('Active set iteration n. %i \n',itAS)
               % Update dof based on current active set
               obj.updateDoFMap;

               % Compute rhs terms (no BCS)
               rhs = computeRhs(obj);

               % assemble global jacobian
               J = computeJacobian2(obj);


               % apply BCs to global system
               [J,rhs] = applyBCFaults(obj,J,rhs,obj.t);


               % compute base rhs norm for relative convergence
               % Norm of unbalanced forces
               rhsNorm = norm(rhs,2);
               rhsNorm0 = rhsNorm;
               itNR = 0;
               k = k+1;
               % update results structure
               obj.results(obj.tStep).itNR(k) = itNR;
               obj.results(obj.tStep).itAS(k) = itAS;
               obj.results(obj.tStep).rhsNorm(k) = 1;
               % NR printing
               if obj.simParameters.verbosity > 1
                  fprintf('Iter     ||rhs||\n');
               end
               % Relative residual
               tolWeigh = obj.simParameters.relTol*rhsNorm;


               if obj.simParameters.verbosity > 1
                  fprintf('0     %e\n',rhsNorm);
               end
               % Active set loop
               while ((rhsNorm > tolWeigh) && (itNR < obj.simParameters.itMaxNR) ...
                     && (rhsNorm > absTol)) || itNR == 0
                  itNR = itNR + 1;
                  %
                  k = k+1; 
                  du = J\(-rhs);

                  % update solution fields (disp, multipliers)
                  updateStateFaults(obj,du);

                  % update nodal gaps
                  dofStick = get_dof(obj.activeSet.curr.stick);
                  computeNodalGap(obj.mortar,obj.state,obj.dofMap,obj.currMultipliers,dofStick);

                  % get information on slipping nodes
                  % if ~isempty(obj.activeSet.curr.slip)
                  %    % get gap from extreme nodes of the fault
                  %    nm = [6;7]; % boundary master nodes at top
                  %    ns = [5;8]; % boundary slave nodes at top
                  %    dum = du(obj.mortar.getContactDofs('master',nm));
                  %    um = obj.state(1).dispCurr(get_dof(nm));
                  %    dus = du(obj.mortar.getContactDofs('slave',ns));
                  %    us = obj.state(2).dispCurr(get_dof(ns));
                  %    fprintf('y = 0 \n')
                  %    fprintf('master du = %.3e %.3e %.3e, slave du = %.3e %.3e %.3e \n',...
                  %       dum(1:3)', dus(1:3)');
                  %    fprintf('master disp = %.3e %.3e %.3e, slave disp = %.3e %.3e %.3e \n',um(1:3)',us(1:3)')
                  %    fprintf('y = 10 \n')
                  %    fprintf('master du = %.3e %.3e %.3e, slave du = %.3e %.3e %.3e \n',...
                  %       dum(4:6)', dus(4:6)');
                  %    fprintf('master disp = %.3e %.3e %.3e, slave disp = %.3e %.3e %.3e \n',um(4:6)',us(4:6)')
                  % end

                  % compute mechanics matrices (linear)
                  for i = 1:obj.nDom
                     getPoro(obj.models(i).Discretizer).computeMat(obj.state(i),obj.dt);
                  end

                  % recompute rhs
                  rhs = computeRhs(obj);

                  % assemble jacobian
                  J = computeJacobian2(obj);

                  % Apply BCs to global system
                  [J,rhs] = applyBCFaults(obj,J,rhs,obj.t);

                  % compute Rhs norm
                  %rhsNorm = norm(rhs(1:3*(obj.mortar.totNodMaster+obj.mortar.totNodSlave)),2);
                  rhsNorm = norm(rhs,2);
                  if obj.simParameters.verbosity > 1
                     fprintf('%d     %e\n',itNR,rhsNorm);
                  end
                  % update results structure
                  obj.results(obj.tStep).itNR(k) = itNR;
                  obj.results(obj.tStep).itAS(k) = itAS;
                  obj.results(obj.tStep).rhsNorm(k) = rhsNorm/rhsNorm0; 

               end % end newton
               %
               % Check NR convergence
               flConv = (rhsNorm < tolWeigh || rhsNorm < absTol);
               if flConv % Convergence
                  for i = 1:obj.nDom
                     obj.state(i).t = obj.t;
                     obj.state(i).advanceState();
                  end
                  % store previous active set
                  obj.activeSet.prev = obj.activeSet.curr;
                  flagActiveSet = updateActiveSet(obj);
                  itAS = itAS+1;
                  if itAS > obj.maxActiveSetIters
                     fprintf('ACTIVE SET ITERATION NOT SUFFICIENT \n')
                     flagActiveSet = false;
                  end
               end
               % Print converged active set solution
               if (flConv) && (~flagActiveSet)
                  if obj.t > obj.simParameters.tMax   % For Steady State
                     for i = 1:obj.nDom
                        printState(obj.models(i).OutState,obj.state(i));
                     end
                     printID = printFault(obj,tListFault,vtkFault,printID);
                  else
                     for i = 1:obj.nDom
                        printState(obj.models(i).OutState,obj.stateOld(i),obj.state(i));
                     end
                     printID = printFault(obj,tListFault,vtkFault,printID,'transient');
                  end
                  sn = obj.currMultipliers(3*id1-2);
                  tauNorm = sqrt(obj.currMultipliers(3*id1).^2+obj.currMultipliers(3*id1-1).^2);
                  uz = obj.mortar.gap(3*id1);
                  obj.results(obj.tStep).s_n = sn;
                  obj.results(obj.tStep).tauNorm = tauNorm;
                  obj.results(obj.tStep).gap = uz;
               end
               %
               manageNextTimeStep(obj,flConv,flagActiveSet);
            end
            %
         end
         %
         vtkFault.finalize();
      end

      function setNonLinearSolverFaults(obj,mG)
         obj.models = mG.model;
         % add lagrange multiplier field to state structure
         obj.currMultipliers = zeros(3*obj.mortar.nMult,1);
         obj.prevMultipliers = obj.currMultipliers;
         obj.mortar.gap = zeros(3*obj.mortar.nS,1);  % NO INITIAL GAP
         % state struct easier to access (we only work with a curr. state)
         obj.state = [obj.models(1).State;obj.models(2).State];
         obj.stateOld = [copy(obj.models(1).State);copy(obj.models(2).State)];
         setInitialStress(obj,'x',-1);         % initialize stress (not flexible at all)
         % initialize active set structure
         obj.activeSet = struct('prev',[],'curr',[]);
         obj.activeSet.prev = struct('stick',(1:obj.mortar.nMult)','slip',[],'open',[]);
         obj.activeSet.curr = obj.activeSet.prev;
         % define dirichlet nodes (always stick)
         obj.dirichNodes = find(obj.mortar.meshGlue.interfaces.mortar.intSlave.coordinates(:,3) == 0);
         % build dof map with initial active set status
         buildDoFMap(obj);
      end


      function manageNextTimeStep(obj,flConv,flAS)
         % flConv -> newton raphson converged
         % flActSet -> active set changed
         if ~flConv
            for i = 1:obj.nDom
               transferState(obj.stateOld(i),obj.state(i));
            end
            obj.t = obj.t - obj.dt;
            obj.tStep = obj.tStep - 1;
            obj.dt = obj.dt/obj.simParameters.divFac;  % Time increment chop
            if obj.dt < obj.simParameters.dtMin
               error('Minimum time step reached')
            elseif obj.simParameters.verbosity > 0
               fprintf('\n %s \n','BACKSTEP');
            end
            % go back to previous active-set
            obj.activeSet = obj.convActiveSet;
            % recompute the nodal gap
            computeNodalGap(obj.mortar,obj.state,obj.dofMap,obj.currMultipliers,get_dof(obj.activeSet.curr.stick));
            return
         elseif flConv && ~flAS
            % update time step
            obj.dt = min([obj.dt*obj.simParameters.multFac,obj.simParameters.dtMax]);
            obj.dt = max([obj.dt obj.simParameters.dtMin]);
            for i = 1:obj.nDom
               transferState(obj.state(i),obj.stateOld(i));
            end
            obj.prevMultipliers = obj.currMultipliers;
            % save active set for possible backstep
            obj.convActiveSet = obj.activeSet;
            %
            if ((obj.t + obj.dt) > obj.simParameters.tMax)
               obj.dt = obj.simParameters.tMax - obj.t;
            end
         end
      end

      function updateStateFaults(obj,du)
         % update solution vectors and derived quantities (strains)
         obj.state(obj.mortar.tagMaster).dispCurr = obj.state(obj.mortar.tagMaster).dispCurr + ...
            du(obj.dofMap.master);
         obj.state(obj.mortar.tagSlave).dispCurr = obj.state(obj.mortar.tagSlave).dispCurr + ...
            du(obj.dofMap.slave);
         obj.currMultipliers = obj.currMultipliers + du(obj.dofMap.slave(end)+1:end);
         %
         % Update stress
         for i = 1:obj.nDom
            l = 0;
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


      function rhsVal = computeRhsMechanics(obj)
         % Compute rhs of individual domains
         for i = 1:obj.nDom
            % Update tmpState
            getPoro(obj.models(i).Discretizer).computeRhs(obj.state(i));
         end
         % residual standard internal forces
         rhsMaster = getField(obj.models(obj.mortar.tagMaster).Discretizer,'Poro').rhs;
         rhsSlave = getField(obj.models(obj.mortar.tagSlave).Discretizer,'Poro').rhs;
         rhsVal = [rhsMaster; rhsSlave];
      end

      function buildDoFMap(obj)
         obj.dofMap = struct('master',[],...
            'slave',[],'intMaster',[],'intSlave',[],'stick',[],'slip',[],...
            'open',[],'nodSlave',[],'nodMaster',[]);
         obj.dofMap.master = getContactDofs(obj,'master');
         obj.dofMap.slave = getContactDofs(obj,'slave');
         obj.dofMap.intMaster = getContactDofs(obj,'master',obj.mortar.idMaster');
         obj.dofMap.intSlave = getContactDofs(obj,'slave',obj.mortar.idSlave');
         obj.dofMap.stick = getContactDofs(obj,'lag',obj.activeSet.curr.stick);
         obj.dofMap.slip = getContactDofs(obj,'lag',obj.activeSet.curr.slip);
         obj.dofMap.open = getContactDofs(obj,'lag',obj.activeSet.curr.open);
         obj.dofMap.nodSlave = get_dof(obj.mortar.meshGlue.interfaces.mortar.nodesSlave);
         obj.dofMap.nodMaster = get_dof(obj.mortar.meshGlue.interfaces.mortar.nodesMaster);
      end

      function updateDoFMap(obj)
         obj.dofMap.stick = getContactDofs(obj,'lag',obj.activeSet.curr.stick);
         obj.dofMap.slip = getContactDofs(obj,'lag',obj.activeSet.curr.slip);
         obj.dofMap.open = getContactDofs(obj,'lag',obj.activeSet.curr.open);
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
                  dofIn = 1:obj.mortar.totNodMaster;
               end
               dofs = get_dof(dofIn);
            case 'slave'
               if isempty(dofIn)
                  dofIn = 1:obj.mortar.totNodSlave;
               end
               dofs = 3*obj.mortar.totNodMaster+get_dof(dofIn);
            case 'lag'
               if isempty(dofIn)
                  dofs = [];
                  return
               end
               dofs = 3*(obj.mortar.totNodMaster+obj.mortar.totNodSlave)+get_dof(dofIn);
            otherwise
               error(['Invalide tag string for getContactDofs method \n:' ...
                  'valid inputs are: master, slave, lag']);
         end
      end

      function rhs = computeRhs(obj)
         % global rhs computation method
         % rhs contribution of BCs is considered in another method!
         rhs = zeros(obj.mortar.nDoF,1);
         % Mechanics internal forces
         rhsMech = computeRhsMechanics(obj);
         rhs([obj.dofMap.master;obj.dofMap.slave]) = rhs([obj.dofMap.master;obj.dofMap.slave])+ rhsMech;
         % Mesh tying rhs
         rhsMeshTying = computeRhsMeshTying(obj);
         rhs([obj.dofMap.intMaster;obj.dofMap.intSlave]) = ...
            rhs([obj.dofMap.intMaster;obj.dofMap.intSlave]) + rhsMeshTying;
         % Stick nodes rhs
         rshStick = computeRhsStick(obj);
         rhs(obj.dofMap.stick) = rhs(obj.dofMap.stick) + rshStick;
         if strcmp(obj.multType,'P0')
            dof = [obj.dofMap.stick; obj.dofMap.slip(1:3:end)];
            rhsStab = computeRhsStabilization(obj);
            rhs(dof) = rhs(dof) + rhsStab;
         end
         % Slip nodes rhs
         rhsSlip = computeRhsSlip(obj);
         rhs(obj.dofMap.slip) = rhs(obj.dofMap.slip) + rhsSlip;
         % Open nodes rhs
         rhsOpen = computeRhsOpen(obj);
         rhs(obj.dofMap.open) = rhs(obj.dofMap.open) + rhsOpen;
         %
         fprintf('Rhs mechanics: %.3e \n',norm(rhsMech));
         fprintf('Rhs mesh tying: %.3e \n',norm(rhsMeshTying))
         fprintf('Rhs stick: %.3e \n',norm(rshStick))
         fprintf('Rhs slip: %.3e \n',norm(rhsSlip))
         fprintf('Rhs open: %.3e \n',norm(rhsOpen))
      end

      function J = computeJacobian(obj)
         % dof ordering: master - slave - lagrange
         % interface dofs are not separated by inner dofs numbering
         i = []; j = []; Jvec = [];
         J = sparse(obj.mortar.nDoF,obj.mortar.nDoF);
         lagDof = getContactDofs(obj,'lag',1:obj.mortar.nMult); % complete set of multipliers
         % mechanics block
         Km = getPoro(obj.models(obj.mortar.tagMaster).Discretizer).K;
         Ks = getPoro(obj.models(obj.mortar.tagSlave).Discretizer).K;
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.master,obj.dofMap.master,Km);
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slave,obj.dofMap.slave,Ks);
         clear Km; clear Ks;
         % mesh tying blocks
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.intMaster,lagDof,-obj.mortar.Mg');
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.intSlave,lagDof,obj.mortar.Dg');
         % stick blocks
         stickDofs = get_dof(obj.activeSet.curr.stick);
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.stick,obj.dofMap.intMaster,-obj.mortar.Mg(stickDofs,:));
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.stick,obj.dofMap.intSlave,obj.mortar.Dg(stickDofs,:));
         % [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.stick,obj.dofMap.intMaster,-obj.mortar.Mn(stickDofs,:));
         % [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.stick,obj.dofMap.intSlave,obj.mortar.Dn(stickDofs,:));
         % [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.stick,obj.dofMap.intMaster,-obj.mortar.Mt(stickDofs,:));
         % [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.stick,obj.dofMap.intSlave,obj.mortar.Dt(stickDofs,:));
         % slip blocks
         if any(obj.activeSet.curr.slip)
            slipDofs = get_dof(obj.activeSet.curr.slip);
            % contribution to normal component
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.intMaster,-obj.mortar.Mn(slipDofs,:));
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.intSlave,obj.mortar.Dn(slipDofs,:));
            % consistency matrices
            [TD,TM,N] = obj.mortar.computeConsistencyMatrices(obj.activeSet,obj.currMultipliers);
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.intMaster,TM(slipDofs,:));
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.intSlave,-TD(slipDofs,:));
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.slip,-N(slipDofs,slipDofs));
            tComp = repmat([false;true;true],numel(obj.activeSet.curr.slip),1);
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.slip,...
               tComp'.*obj.mortar.L(slipDofs,slipDofs).*tComp);

            % % contribution to tangential component
            % T = computeDtDgt(obj.mortar,obj.activeSet,obj.currMultipliers); % non linear contribution
            % [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.intMaster,T*obj.mortar.Mt(slipDofs,:));
            % [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.intSlave,-T*obj.mortar.Dt(slipDofs,:));
            % %N = computeDtDtn(obj.mortar,obj.activeSet,obj.currMultipliers); % non linear contribution
            % %select only tangential components of mass matrix
            % tComp = repmat([false;true;true],numel(obj.activeSet.curr.slip),1);
            % %provisional matrices consistent with Jha formulation
            % % [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.slip,...
            % %    -N*obj.mortar.L(slipDofs,slipDofs));
            % [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.slip,...
            %    tComp'.*obj.mortar.L(slipDofs,slipDofs).*tComp);
         end
         % open blocks
         if any(obj.activeSet.curr.open)
            openDofs = get_dof(obj.activeSet.curr.open);
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.open,obj.dofMap.open,...
               obj.mortar.L(openDofs,openDofs));
         end
         J = sparse(i,j,Jvec,obj.mortar.nDoF,obj.mortar.nDoF);
      end

      function J = computeJacobian2(obj)
         % dof ordering: master - slave - lagrange
         % interface dofs are not separated by inner dofs numbering
         J = sparse(obj.mortar.nDoF,obj.mortar.nDoF);
         lagDof = getContactDofs(obj,'lag',1:obj.mortar.nMult); % complete set of multipliers
         % mechanics block
         Km = getPoro(obj.models(obj.mortar.tagMaster).Discretizer).K;
         Ks = getPoro(obj.models(obj.mortar.tagSlave).Discretizer).K;
         J(obj.dofMap.master,obj.dofMap.master) = Km;
         J(obj.dofMap.slave,obj.dofMap.slave) = Ks;
         clear Km; clear Ks;
         % mesh tying blocks
         J(obj.dofMap.intMaster,lagDof) = -obj.mortar.Mg';
         J(obj.dofMap.intSlave,lagDof) = obj.mortar.Dg';
         % stick blocks
         stickDofs = get_dof(obj.activeSet.curr.stick);
         J(obj.dofMap.stick,obj.dofMap.intMaster) = J(obj.dofMap.stick,obj.dofMap.intMaster) - obj.mortar.Mn(stickDofs,:);
         J(obj.dofMap.stick,obj.dofMap.intSlave) = J(obj.dofMap.stick,obj.dofMap.intSlave) + obj.mortar.Dn(stickDofs,:);
         J(obj.dofMap.stick,obj.dofMap.intMaster) = J(obj.dofMap.stick,obj.dofMap.intMaster) - obj.mortar.Mt(stickDofs,:);
         J(obj.dofMap.stick,obj.dofMap.intSlave) = J(obj.dofMap.stick,obj.dofMap.intSlave) + obj.mortar.Dt(stickDofs,:);
         if strcmp(obj.multType,'P0')
            % stabilization contribution to stick dofs
            J(obj.dofMap.stick,obj.dofMap.stick) = J(obj.dofMap.stick,obj.dofMap.stick) - obj.mortar.stabMat(stickDofs,stickDofs);
         end
         % slip blocks
         if any(obj.activeSet.curr.slip)
            % split normal and tangential components
            slipDofs = get_dof(obj.activeSet.curr.slip);
            if strcmp(obj.multType,'P0')
               slipH = slipDofs(1:3:end);
               slipJ = obj.dofMap.slip;
               slipJ = slipJ(1:3:end);
               % stabilization contribution to slip normal components and
               % adjacent components
               J(obj.dofMap.stick,slipJ) = J(obj.dofMap.stick,slipJ) - obj.mortar.stabMat(stickDofs,slipH);
               J(slipJ,obj.dofMap.stick) = J(slipJ,obj.dofMap.stick) - obj.mortar.stabMat(slipH,stickDofs);
               J(slipJ,slipJ) = J(slipJ,slipJ) - obj.mortar.stabMat(slipH,slipH);
            end
            % contribution to normal component
            J(obj.dofMap.slip,obj.dofMap.intMaster) = -obj.mortar.Mn(slipDofs,:);
            J(obj.dofMap.slip,obj.dofMap.intSlave) = obj.mortar.Dn(slipDofs,:);
            % consistency matrices
            [TD,TM,N] = obj.mortar.computeConsistencyMatrices(obj.activeSet,obj.currMultipliers,obj.state,obj.stateOld,obj.dofMap);
            J(obj.dofMap.slip,obj.dofMap.intMaster) = J(obj.dofMap.slip,obj.dofMap.intMaster)+TM(slipDofs,:);
            J(obj.dofMap.slip,obj.dofMap.intSlave) = J(obj.dofMap.slip,obj.dofMap.intSlave)-TD(slipDofs,:);
            J(obj.dofMap.slip,obj.dofMap.slip) = J(obj.dofMap.slip,obj.dofMap.slip)-N(slipDofs,slipDofs);
            tComp = repmat([false;true;true],numel(obj.activeSet.curr.slip),1);
            J(obj.dofMap.slip,obj.dofMap.slip) = J(obj.dofMap.slip,obj.dofMap.slip)+tComp'.*obj.mortar.L(slipDofs,slipDofs).*tComp;
         end
         % open blocks
         if any(obj.activeSet.curr.open)
            openDofs = get_dof(obj.activeSet.curr.open);
            J(obj.dofMap.open,obj.dofMap.open) = obj.mortar.L(openDofs,openDofs);
         end
      end



      function rhsVal = computeRhsMeshTying(obj)
         % subtract initial stress state
         rhsMaster = -obj.mortar.Mg'*(obj.currMultipliers-obj.iniMult);
         rhsSlave = obj.mortar.Dg'*(obj.currMultipliers-obj.iniMult);
         rhsVal = [rhsMaster;rhsSlave];
      end

      function isChanged = updateActiveSet(obj)
         % update and compare active set
         % dirichlet nodes are not subject to contact check and remain
         % stick
         %idDir = ismember(obj.activeSet.curr.stick,obj.dirichNodes);
         % return true if some changes are made
         % return false if the active set stays the same
         multMat = (reshape(obj.currMultipliers,3,[]))';
         % module of limit traction
         tLim = obj.mortar.coes - multMat(:,1)*tan(deg2rad(obj.mortar.phi));
         tauNorm = sqrt(multMat(:,2).^2+multMat(:,3).^2);
         % stick mode to slip mode
         stick2slip = tauNorm(obj.activeSet.curr.stick) > (1+obj.mortar.tolTang)*tLim(obj.activeSet.curr.stick);
         obj.activeSet.curr.slip = unique([obj.activeSet.curr.slip; ...
             obj.activeSet.curr.stick(stick2slip)]);
         % stick mode to open mode 
         stick2open = multMat(obj.activeSet.curr.stick,1) > obj.mortar.tolNormal;
         % stick2open(idDir) = false;
         obj.activeSet.curr.open = unique([obj.activeSet.curr.open; ...
             obj.activeSet.curr.stick(stick2open)]);
         % open mode to stick node
         open2stick = obj.mortar.g_N(obj.activeSet.curr.open) < obj.mortar.tolGap;
         id = all([stick2open stick2slip],2);
         stick2open(id) = false;
         obj.activeSet.curr.stick = unique([obj.activeSet.curr.stick(~any([stick2open stick2slip],2));...
             obj.activeSet.curr.open(open2stick)]);
         obj.activeSet.curr.open = obj.activeSet.curr.open(~open2stick);
         % compare current and previous active set
         isChanged = compareSet(obj.activeSet.curr.open,obj.activeSet.prev.open);
         if isChanged
            return
         else
            isChanged = compareSet(obj.activeSet.curr.stick,obj.activeSet.prev.stick);
         end
         if isChanged
            return
         else
            isChanged = compareSet(obj.activeSet.curr.slip,obj.activeSet.prev.slip);
         end
      end

      function initActiveSet(obj)
          % move slip nodes to stick condition
          obj.activeSet.curr.stick = unique([obj.activeSet.curr.stick;...
              obj.activeSet.curr.slip]);
          obj.activeSet.curr.slip = [];
      end

      function rhsStick = computeRhsStick(obj)
         usCurr = obj.state(obj.mortar.tagSlave).dispCurr(obj.dofMap.nodSlave);
         usConv = obj.state(obj.mortar.tagSlave).dispConv(obj.dofMap.nodSlave);
         umCurr = obj.state(obj.mortar.tagMaster).dispCurr(obj.dofMap.nodMaster);
         umConv = obj.state(obj.mortar.tagMaster).dispConv(obj.dofMap.nodMaster);
         rhsStick = obj.mortar.Dn*usCurr - obj.mortar.Mn*umCurr + ...
            obj.mortar.Dt*(usCurr-usConv) - obj.mortar.Mt*(umCurr-umConv);
         %rhsStick = obj.mortar.Dg*usCurr - obj.mortar.Mg*umCurr;
         rhsStick = rhsStick(get_dof(obj.activeSet.curr.stick));
         %rhsStick = 0;
      end

      function rhsSlip = computeRhsSlip(obj)
         if ~isempty(obj.activeSet.curr.slip)
            dofSlip = get_dof(obj.activeSet.curr.slip);
            tracLim = computeLimitTraction(obj.mortar,obj.activeSet,obj.currMultipliers); % limit traction in active nodes
            usCurr = obj.state(obj.mortar.tagSlave).dispCurr(obj.dofMap.nodSlave);
            umCurr = obj.state(obj.mortar.tagMaster).dispCurr(obj.dofMap.nodMaster);
            % if I have an available non-zero slip, correct the current
            % traction direction according to the direction of the slip
            % get sign of the slip
            %slipDir = obj.mortar.g_T(dofSlip)./abs(obj.mortar.g_T(dofSlip));
            %slipDir(isnan(slipDir)) = 0; % normal direction does not count
            %slipCorr = abs(obj.mortar.g_T(dofSlip)) > 1e-10;  
            rhsSlip1 = obj.mortar.Dn*usCurr - obj.mortar.Mn*umCurr; % normal components
            t_T = repmat([false;true;true],numel(obj.activeSet.curr.slip),1).*(obj.currMultipliers(dofSlip));
            %t_T(slipCorr) = slipDir(slipCorr).*abs(t_T(slipCorr));
            %tracError = norm(t_T - tracLim);
            rhsSlip2 = obj.mortar.L(dofSlip,dofSlip)*t_T; % tangential components
            %rhsTest =  obj.mortar.L(dofSlip,dofSlip)*tracLim;
            %rhsSlip3 = computeRhsLimitTraction(obj.mortar,obj.currMultipliers,obj.activeSet);
            rhsSlip4 = computeRhsLimitTraction2(obj.mortar,...
               obj.currMultipliers,obj.state,obj.stateOld,obj.activeSet,obj.dofMap);
            % if obj.t>5.5
            %    rhsSlip2 = -rhsSlip2;
            %    if ~any(usCurr)
            %       rhsSlip4 = -rhsSlip4;
            %    end
            % end
            rhsSlip = rhsSlip1(dofSlip) + rhsSlip2 - rhsSlip4(dofSlip);
         else
            rhsSlip = [];
         end
      end

      function rhsOpen = computeRhsOpen(obj)
         if ~isempty(obj.activeSet.curr.open)
            rhsOpen = obj.mortar.L*obj.currMultipliers;
            rhsOpen = rhsOpen(get_dof(obj.activeSet.curr.open));
         else
            rhsOpen = [];
         end
      end

      function rhsStab = computeRhsStabilization(obj)
         % compute the stabilization contribution to the rhs
         if ~strcmp(obj.multType,'P0')
            dof = get_dof(obj.activeSet.curr.stick);
         else
            dofS = get_dof(obj.activeSet.curr.stick);
            dofN = get_dof(obj.activeSet.curr.slip);
            dofN = dofN(1:3:end);
            dof = [dofS;dofN];
         end
         rhsStab = -obj.mortar.stabMat(dof,dof)*(obj.currMultipliers(dof)-obj.iniMult(dof));
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
         v_loc = obj.mortar.R*repmat(v',obj.mortar.nMult,1);
         % make sure contact traction is negative (compression)
         obj.currMultipliers(:) = -abs(v_loc);
         obj.iniMult = obj.currMultipliers;
         obj.state(1).conv.stress = obj.state(1).iniStress;
         %obj.state(1).curr.stress = obj.state(1).iniStress;
         obj.state(2).conv.stress = obj.state(2).iniStress;
         %obj.state(2).curr.stress = obj.state(2).iniStress;
      end
   end

   methods (Access=private)
      function tID = printFault(obj,tList,vtk,tID,var)
         if nargin == 4
            time = obj.stateOld(1).t;
            outVar = obj.prevMultipliers;
            data = repmat(struct('name',[],'data',[]),3,1);
            name = ["sigma_n","tau_1","tau_2"];
            for i = 1:size(outVar,2)
               data(i).name = convertStringsToChars(name(i));
               data(i).data = outVar(:,i);
            end
            vtk.writeVTKFile(time, [], [], data, []);
         elseif nargin == 5
            while tID <= length(tList) &&(tList(tID)<= obj.state(1).t)
               % assert(tList(tID) > obj.stateOld(1).t, ...
               %    'Print time %f out of range (%f - %f)',timeList(tID), ...
               %    obj.stateOld(1).t,obj.state(1).t);
               assert(obj.state(1).t - obj.stateOld(1).t > eps('double'),'Dt too small for printing purposes');
               %
               % Linear interpolation
               fac = (tList(tID) - obj.stateOld(1).t)/(obj.state(1).t - obj.stateOld(1).t);
               time = tList(tID);
               outVar = zeros(obj.mortar.nMult,6);
               interpMult= fac*obj.currMultipliers + (1-fac)*obj.prevMultipliers;
               outVar(:,1:3) = (reshape(interpMult,3,[]))';
               outVar(:,4) = sqrt(outVar(:,2).^2 + outVar(:,3).^2);
               % gap values
               outVar(:,5) = obj.mortar.g_N;
               gt = (reshape(obj.mortar.g_T,3,[]))';
               outVar(:,6) = sqrt(gt(:,2).^2+gt(:,3).^2);
               data = repmat(struct('name',[],'data',[]),3,1);
               name = ["sigma_n","tau_1","tau_2","tau_norm","g_N","||g_T||"];
               for i = 1:size(outVar,2)
                  data(i).name = convertStringsToChars(name(i));
                  data(i).data = outVar(:,i);
               end
               if strcmp(obj.multType,'P0')
                  vtk.writeVTKFile(time, [], [], [], data);
               else
                  vtk.writeVTKFile(time, [], [], data, []);
               end
               tID = tID + 1;
            end
         end
      end
   end
end
