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
         vtkFault = VTKOutput(obj.mortar.meshGlue.interfaces.mortar.intSlave,'Fault');
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
            obj.tStep = obj.tStep + 1;
            obj.t = obj.t + obj.dt;
            % Apply Dirichlet value to individual domain solutions
            % Apply the Dirichlet condition value to the solution vector
            applyDirFault(obj);

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

               % Newton Raphson loop
               while ((rhsNorm > tolWeigh) && (itNR < obj.simParameters.itMaxNR) ...
                     && (rhsNorm > absTol)) || itNR == 0
                  itNR = itNR + 1;
                  %
                  du = J\(-rhs);

                  % update solution fields (disp, multipliers)
                  updateStateFaults(obj,du);

                  % printe state (even if not converged)
                  % for i = 1:obj.nDom
                  %     printState(obj.models(i).OutState,obj.state(i));
                  %     obj.models(i).OutState.finalize();
                  %     plotFunction(obj.mortar.meshGlue.interfaces.mortar.intSlave,'Fault',(reshape(obj.currMultipliers,3,[]))',...
                  %         ["sigma_n","tau_1","tau_2"]);
                  % end

                  % update nodal gaps
                  computeNodalGap(obj.mortar,obj.state,obj.dofMap);

                  % recompute rhs
                  rhs = computeRhs(obj);

                  % assemble jacobian
                  J = computeJacobian(obj);

                  % Apply BCs to global system
                  [J,rhs] = applyBCFaults(obj,J,rhs,obj.t);

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
                     obj.state(i).t = obj.t;
                     % Print the solution, if needed
                     obj.state(i).advanceState();
                  end
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
                  % Update active set
                  flagActiveSet = updateActiveSet(obj);
                  obj.activeSet.prev = obj.activeSet.curr;
               end
            end
            %
            % Manage next time step
            manageNextTimeStep(obj,flConv);
         end
         vtkFault.finalize();
         %
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

      function manageNextTimeStep(obj,flConv)
         if ~flConv   % Perform backstep
            for i = 1:obj.nDom
               transferState(obj.stateOld(i),obj.state(i));
            end
            obj.prevMultipliers = obj.currMultipliers;
            obj.t = obj.t - obj.dt;
            obj.tStep = obj.tStep - 1;
            obj.dt = obj.dt/obj.simParameters.divFac;
            obj.dt = obj.dt/obj.simParameters.divFac;  % Time increment chop
            if obj.dt < obj.simParameters.dtMin
               error('Minimum time step reached')
            elseif obj.simParameters.verbosity > 0
               fprintf('\n %s \n','BACKSTEP');
            end
         end

         if flConv % Go on if converged
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
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.stick,obj.dofMap.intMaster,-obj.mortar.Mn);
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.stick,obj.dofMap.intSlave,obj.mortar.Dn);
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.stick,obj.dofMap.intMaster,-obj.mortar.Mt);
         [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.stick,obj.dofMap.intSlave,obj.mortar.Dt);
         % slip blocks
         if any(obj.activeSet.curr.slip)
            slipDofs = get_dof(obj.activeSet.curr.slip);
            % contribution to normal component
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.intMaster,-obj.mortar.Mn(slipDofs,:));
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.intSlave,obj.mortar.Dn(slipDofs,:));
            % contribution to tangential component
            T = computeDtDgt(obj); % non linear contribution
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.intMaster,-T*obj.mortar.Mt(slipDofs,:));
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.intSlave,T*obj.mortar.Dt(slipDofs,:));
            N = computeDtDtn(obj); % non linear contribution
            % select only tangential components of mass matrix
            tComp = repmat([false;true;true],numel(obj.activeSet.curr.slip),1);
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.slip,...
               N*obj.mortar.L(slipDofs,slipDofs));
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.slip,obj.dofMap.slip,...
               tComp'.*obj.mortar.L(slipDofs,slipDofs).*tComp);
         end
         % open blocks
         if any(obj.activeSet.curr.open)
            openDofs = get_dof(obj.activeSet.curr.open);
            [i,j,Jvec] = addBlockJ(i,j,Jvec,obj.dofMap.open,obj.dofMap.open,...
               obj.mortar.L(openDofs,openDofs));
         end
         J = sparse(i,j,Jvec,obj.mortar.nDoF,obj.mortar.nDoF);
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
         rhsStick = rhsStick(get_dof(obj.activeSet.curr.stick));
      end

      function rhsSlip = computeRhsSlip(obj)
         if ~isempty(obj.activeSet.curr.slip)
            dofSlip = get_dof(obj.activeSet.curr.slip);
            tracLim = computeLimitTraction(obj); % limit traction in active nodes
            usCurr = obj.state(obj.mortar.tagSlave).dispCurr(obj.dofMap.nodSlave);
            umCurr = obj.state(obj.mortar.tagMaster).dispConv(obj.dofMap.nodMaster);
            rhsSlip1 = obj.mortar.Dn*usCurr - obj.mortar.Mn*umCurr; % normal components
            rhsSlip2 = obj.mortar.L(dofSlip,dofSlip)*(obj.currMultipliers(dofSlip)-tracLim); % tangential components
            rhsSlip = rhsSlip1(dofSlip) + rhsSlip2;
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
         obj.state(2).conv.stress = obj.state(2).iniStress;
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
            if tID <= length(tList)
               while (tList(tID)<= obj.state(1).t)
                  % assert(tList(tID) > obj.stateOld(1).t, ...
                  %    'Print time %f out of range (%f - %f)',timeList(tID), ...
                  %    obj.stateOld(1).t,obj.state(1).t);
                  assert(obj.state(1).t - obj.stateOld(1).t > eps('double'),'Dt too small for printing purposes');
                  %
                  % Linear interpolation
                  fac = (tList(tID) - obj.stateOld(1).t)/(obj.state(1).t - obj.stateOld(1).t);
                  time = tList(tID);
                  outVar = fac*obj.currMultipliers + (1-fac)*obj.prevMultipliers;
                  outVar = (reshape(outVar,3,[]))';
                  pointData = repmat(struct('name',[],'data',[]),3,1);
                  name = ["sigma_n","tau_1","tau_2"];
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
end
