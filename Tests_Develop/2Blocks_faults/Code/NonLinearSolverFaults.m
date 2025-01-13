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
      maxActiveSetIters = 1 % maximum number of active set iterations
      activeSet
      iniMult; % initial multiplier vector
      currMultipliers % store current and previous contact tractions
      prevMultipliers
      dofMap      % group DoF based on master/slave/state
   end

   methods
      function obj = NonLinearSolverFaults(simParam,mG,c,phi)
         obj.simParameters = simParam;
         obj.mortar = mortarFaults(mG,c,phi);
         obj.setNonLinearSolverFaults(mG);
         obj.simulationLoop();
      end

      function simulationLoop(obj)
         % provisional print utilities for faults
         if strcmp(obj.mortar.meshGlue.interfaces.multType,'dual')
            vtkFault = VTKOutput(obj.mortar.meshGlue.interfaces.mortar.intSlave,'Fault_dual');
         elseif strcmp(obj.mortar.meshGlue.interfaces.multType,'standard')
            vtkFault = VTKOutput(obj.mortar.meshGlue.interfaces.mortar.intSlave,'Fault_standard');
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

         % Initialize nodal gap
         computeNodalGap(obj.mortar,obj.state,obj.dofMap);

         % TIME LOOP
         while obj.t < obj.simParameters.tMax
            % Update the simulation time and time step ID
            absTol = obj.simParameters.absTol;

            itAS = 0;
            flagActiveSet = false;

            obj.mortar.gapOld = obj.mortar.gap;

            obj.tStep = obj.tStep + 1;
            obj.t = obj.t + obj.dt;
            % Apply Dirichlet value to individual domain solutions

            if (obj.simParameters.verbosity > 0) && (itAS == 0)
               fprintf('\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,obj.dt);
               fprintf('-----------------------------------------------------------\n');
            end


            % Apply the Dirichlet condition value to the solution vector
            applyDirFault(obj);
            
            while flagActiveSet || itAS==0
               NRrhsNorm = zeros(obj.simParameters.itMaxNR,1);
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
               rhsNorm = norm(rhs(1:3*(obj.mortar.totNodMaster+obj.mortar.totNodSlave)),2);
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
               % Active set loop
               while ((rhsNorm > tolWeigh) && (itNR < obj.simParameters.itMaxNR) ...
                     && (rhsNorm > absTol)) || itNR == 0
                  itNR = itNR + 1;
                  %
                  du = J\(-rhs);

                  % update solution fields (disp, multipliers)
                  updateStateFaults(obj,du);

                  % update nodal gaps
                  computeNodalGap(obj.mortar,obj.state,obj.dofMap);

                  % get information on slipping nodes
                  k = 0;

                  if ~isempty(obj.activeSet.curr.slip)
                     % get gap from extreme nodes of the fault
                     nm = [6;7]; % boundary master nodes at top
                     ns = [5;8]; % boundary slave nodes at top
                     dum = du(obj.mortar.getContactDofs('master',nm));
                     um = obj.state(1).dispCurr(get_dof(nm));
                     dus = du(obj.mortar.getContactDofs('slave',ns));
                     us = obj.state(2).dispCurr(get_dof(ns));
                     fprintf('y = 0 \n')
                     fprintf('master du = %.3e %.3e %.3e, slave du = %.3e %.3e %.3e \n',...
                        dum(1:3)', dus(1:3)');
                     fprintf('master disp = %.3e %.3e %.3e, slave disp = %.3e %.3e %.3e \n',um(1:3)',us(1:3)')
                     fprintf('y = 10 \n')
                     fprintf('master du = %.3e %.3e %.3e, slave du = %.3e %.3e %.3e \n',...
                        dum(4:6)', dus(4:6)');
                     fprintf('master disp = %.3e %.3e %.3e, slave disp = %.3e %.3e %.3e \n',um(4:6)',us(4:6)')
                  end

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

                  NRrhsNorm(itNR) = rhsNorm;

               end % end newton
               %
               % Check NR convergence
               flConv = (rhsNorm < tolWeigh || rhsNorm < absTol);
               if flConv % Convergence
                  for i = 1:obj.nDom
                     obj.state(i).t = obj.t;
                     obj.state(i).advanceState();
                  end
                  % Update active set if NR converged
                  if flagActiveSet
                     figure(1)
                     semilogy(k+1:k+numel(NRrhsNorm),NRrhsNorm,'-o','DisplayName',strcat('Iter_',num2str(itNR)))
                     legend('-DynamicLegend')
                     hold on
                  end
                  k = k + numel(NRrhsNorm);
                  flagActiveSet = updateActiveSet(obj);
                  obj.activeSet.prev = obj.activeSet.curr;
                  itAS = itAS+1;
                  % plot NR convergence
                  NRrhsNorm = NRrhsNorm(NRrhsNorm~=0);
                  if itAS > obj.maxActiveSetIters
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
               end
               %
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
         obj.currMultipliers = zeros(3*obj.mortar.nS,1);
         obj.prevMultipliers = obj.currMultipliers;
         obj.mortar.gap = zeros(3*obj.mortar.nS,1);  % NO INITIAL GAP
         % state struct easier to access (we only work with a curr. state)
         obj.state = [obj.models(1).State;obj.models(2).State];
         obj.stateOld = [copy(obj.models(1).State);copy(obj.models(2).State)];
         setInitialStress(obj,'x',-1);         % initialize stress (not flexible at all)
         % initialize active set structure
         obj.activeSet = struct('prev',[],'curr',[]);
         obj.activeSet.prev = struct('stick',(1:obj.mortar.nS)','slip',[],'open',[]);
         obj.activeSet.curr = obj.activeSet.prev;
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
            return
         elseif flConv && ~flAS
            % update time step
            obj.dt = min([obj.dt*obj.simParameters.multFac,obj.simParameters.dtMax]);
            obj.dt = max([obj.dt obj.simParameters.dtMin]);
            for i = 1:obj.nDom
               transferState(obj.state(i),obj.stateOld(i));
            end
            obj.prevMultipliers = obj.currMultipliers;
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
         lagDof = getContactDofs(obj,'lag',1:obj.mortar.nS); % complete set of multipliers
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
         lagDof = getContactDofs(obj,'lag',1:obj.mortar.nS); % complete set of multipliers
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
         % slip blocks
         if any(obj.activeSet.curr.slip)
            slipDofs = get_dof(obj.activeSet.curr.slip);
            % contribution to normal component
            J(obj.dofMap.slip,obj.dofMap.intMaster) = -obj.mortar.Mn(slipDofs,:);
            J(obj.dofMap.slip,obj.dofMap.intSlave) = obj.mortar.Dn(slipDofs,:);
            % consistency matrices
            [TD,TM,N] = obj.mortar.computeConsistencyMatrices(obj.activeSet,obj.currMultipliers,obj.state,obj.dofMap);
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
         % return true if some changes are made
         % return 0 if the active set stays the same
         list = (1:obj.mortar.nS)';
         multMat = (reshape(obj.currMultipliers,3,[]))';
         tLim = obj.mortar.coes - multMat(:,1)*tan(deg2rad(obj.mortar.phi));
         tauNorm = sqrt(multMat(:,2).^2+multMat(:,3).^2);
         % open mode
         isOpen = any([multMat(:,1) > obj.mortar.tolNormal ...
            obj.mortar.g_N > obj.mortar.tolGap],2);
         obj.activeSet.curr.open = find(isOpen);
         list = list(~isOpen);
         % stick mode
         isStick = tauNorm(list) < tLim(list)*(1+obj.mortar.tolTang);
         obj.activeSet.curr.stick = list(isStick);
         list = list(~isStick);
         % slip mode
         obj.activeSet.curr.slip = list;
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

      function rhsStick = computeRhsStick(obj)
         usCurr = obj.state(obj.mortar.tagSlave).dispCurr(obj.dofMap.nodSlave);
         usConv = obj.state(obj.mortar.tagSlave).dispConv(obj.dofMap.nodSlave);
         umCurr = obj.state(obj.mortar.tagMaster).dispCurr(obj.dofMap.nodMaster);
         umConv = obj.state(obj.mortar.tagMaster).dispConv(obj.dofMap.nodMaster);
         rhsStick = obj.mortar.Dn*usCurr - obj.mortar.Mn*umCurr + ...
            obj.mortar.Dt*(usCurr-usConv) - obj.mortar.Mt*(umCurr-umConv);
         % rhsStick = obj.mortar.Dg*usCurr - obj.mortar.Mg*umCurr;
         rhsStick = rhsStick(get_dof(obj.activeSet.curr.stick));
         %rhsStick = 0;
      end

      function rhsSlip = computeRhsSlip(obj)
         if ~isempty(obj.activeSet.curr.slip)
            dofSlip = get_dof(obj.activeSet.curr.slip);
            tracLim = computeLimitTraction(obj.mortar,obj.activeSet,obj.currMultipliers); % limit traction in active nodes
            usCurr = obj.state(obj.mortar.tagSlave).dispCurr(obj.dofMap.nodSlave);
            umCurr = obj.state(obj.mortar.tagMaster).dispCurr(obj.dofMap.nodMaster);
            rhsSlip1 = obj.mortar.Dn*usCurr - obj.mortar.Mn*umCurr; % normal components
            t_T = repmat([false;true;true],numel(obj.activeSet.curr.slip),1).*obj.currMultipliers(dofSlip);
            tracError = norm(t_T - tracLim);
            rhsSlip2 = obj.mortar.L(dofSlip,dofSlip)*t_T; % tangential components
            rhsTest =  obj.mortar.L(dofSlip,dofSlip)*tracLim;
            rhsSlip3 = computeRhsLimitTraction(obj.mortar,obj.currMultipliers,obj.activeSet);
            rhsSlip4 = computeRhsLimitTraction2(obj.mortar,...
               obj.currMultipliers,obj.state,obj.activeSet,obj.dofMap);
            rhsSlip = rhsSlip1(dofSlip) + rhsSlip2 - rhsSlip4(dofSlip);
         else
            rhsSlip = [];
         end
      end

      function rhsOpen = computeRhsOpen(obj)
         if ~isempty(obj.activeSet.curr.open)
            rhsOpen = obj.mortar.L*obj.currMultipliers;
            rhsOpen = rhsOpen(get_dof(obj.activeSet.curr.slip));
         else
            rhsOpen = [];
         end
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
         v_loc = obj.mortar.R*repmat(v',obj.mortar.nS,1);
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
            pointData = repmat(struct('name',[],'data',[]),3,1);
            name = ["sigma_n","tau_1","tau_2"];
            for i = 1:size(outVar,2)
               pointData(i).name = convertStringsToChars(name(i));
               pointData(i).data = outVar(:,i);
            end
            vtk.writeVTKFile(time, [], [], pointData, []);
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
               outVar = zeros(obj.mortar.nS,6);
               interpMult= fac*obj.currMultipliers + (1-fac)*obj.prevMultipliers;
               outVar(:,1:3) = (reshape(interpMult,3,[]))';
               outVar(:,4) = sqrt(outVar(:,2).^2 + outVar(:,3).^2);
               % gap values
               outVar(:,5) = obj.mortar.g_N;
               gt = (reshape(obj.mortar.g_T,3,[]))';
               outVar(:,6) = sqrt(gt(:,2).^2+gt(:,3).^2);
               pointData = repmat(struct('name',[],'data',[]),3,1);
               name = ["sigma_n","tau_1","tau_2","tau_norm","g_N","||g_T||"];
               for i = 1:size(outVar,2)
                  pointData(i).name = convertStringsToChars(name(i));
                  pointData(i).data = outVar(:,i);
               end
               vtk.writeVTKFile(time, [], [], pointData, []);
               tID = tID + 1;
            end
         end
      end
   end
end
