function computeContactMatricesAndRhs(obj)

% Compute contact matrices and rhs

% the syntax of the comments is inspired by the appendix of
% Franceschini, Andrea, et al. "Algebraically stabilized Lagrange
% multiplier method for frictional contact mechanics with
% hydraulically active fractures." Computer Methods in Applied
% Mechanics and Engineering 368 (2020): 113161.

dofSlave = getDoFManager(obj,MortarSide.slave);
dofMaster = getDoFManager(obj,MortarSide.master);

% extract current and previous displacements
um = obj.domains(1).state.data.displacements;
us = obj.domains(2).state.data.displacements;
um_old = obj.domains(1).stateOld.data.displacements;
us_old = obj.domains(2).stateOld.data.displacements;

% compute displacement increment
%       um_slip = um - um_old;
%       us_slip = us - us_old;

ncomp = dofSlave.getDoFperEnt(obj.physic);

% define matrix assemblers
[asbMu,asbDu,asbMt,asbDt,asbQ] = defineAssemblers(obj);

% define rhs vectors;
rhsUm = zeros(getNumbDoF(dofMaster),1);
rhsUs = zeros(getNumbDoF(dofSlave),1);
rhsT = zeros(getNumbDoF(obj),1);

fldM = dofMaster.getVariableId(obj.coupledVariables);
fldS = dofSlave.getVariableId(obj.coupledVariables);

% anonymous functions for local fem computations
f1 = @(a,b) pagemtimes(a,'ctranspose',b,'none');
f2 = @(a,b,c) pagemtimes(a,'transpose',pagemtimes(b,c),'none');


for iPair = 1:obj.quadrature.numbInterfacePairs

  is = obj.quadrature.interfacePairs(iPair,1);
  im = obj.quadrature.interfacePairs(iPair,2);

  contactState = obj.activeSet.curr(is);

  % retrieve mortar integration data
  xiMaster = obj.quadrature.getMasterGPCoords(iPair);
  xiSlave = obj.quadrature.getSlaveGPCoords(iPair);
  dJw = obj.quadrature.getIntegrationWeights(iPair);

  % define slave related quantities
  nodeSlave = obj.interfMesh.local2glob{2}(obj.interfMesh.msh(2).surfaces(is,:));
  usDof = dofSlave.getLocalDoF(fldS,nodeSlave);
  tDof = getMultiplierDoF(obj,is);
  trac = obj.state.traction(tDof);
  % equilibrium equation and stabilization work with delta traction
  dTrac = trac - obj.iniTraction(tDof);

  % displacement increment to check if this is the first newton
  % iteration
  usDiff = norm(us(usDof) - us_old(usDof));


  nodeMaster = obj.interfMesh.local2glob{1}(obj.interfMesh.msh(1).surfaces(im,:));
  umDof = dofMaster.getLocalDoF(fldM,nodeMaster);

  [Nslave,Nmaster,Nmult] = ...
    getMortarBasisFunctions(obj,im,is,xiMaster,xiSlave);

  % reshape basis function matrices to match number of components
  [Ns,Nm,Nmult] = ...
    obj.reshapeBasisFunctions(ncomp,Nslave,Nmaster,Nmult);

  % rotation of multipliers
  R = getRotationMatrix(obj.interfMesh,is);
  NmultR = pagemtimes(Nmult,R);

  % Reduce the dimension of multiplier basis functions exploiting
  % the local definition of degrees of freedom

  % the normal component of the multipliers basis
  Nmult_n = -Nmult(1,1,:);

  % tangential component of multiplier basis
  Nmult_t = NmultR(:,2:3,:);

  % get normal at the gauss points (for warped facets...?)
  normalNodes = obj.interfMesh.getNodalNormal(is);
  normal = pagemtimes(Ns,normalNodes);

  % operator selecting only tangential components of the
  % displacements
  T = eye(3) - pagemtimes(normal,'none',normal,'transpose');
  % normal and tangential component of displacement basis functions
  Nm_n = pagemtimes(normal,'transpose',Nm,'none');
  Nm_t = pagemtimes(T,Nm);

  Ns_n = pagemtimes(normal,'transpose',Ns,'none');
  Ns_t = pagemtimes(T,Ns);


  % tangential gap increment (quasi-static coulomb low)
  gt_curr = pagemtimes(Ns_t,us(usDof)) ...
    - pagemtimes(Nm_t,um(umDof));
  gt_old =  pagemtimes(Ns_t,us_old(usDof)) ...
    - pagemtimes(Nm_t,um_old(umDof));
  dgt = gt_curr - gt_old;


  % normal gap (scalar)
  g_n = pagemtimes(Ns_n,us(usDof)) ...
    - pagemtimes(Nm_n,um(umDof));


  % A_us
  Aun_m =  MortarQuadrature.integrate(f1,Nm_n,Nmult_n,dJw);
  Aun_s =  MortarQuadrature.integrate(f1,Ns_n,Nmult_n,dJw);
  asbMu.localAssembly(umDof,tDof(1),-Aun_m);
  asbDu.localAssembly(usDof,tDof(1),Aun_s);

  Aut_m =  MortarQuadrature.integrate(f1,Nm_t,Nmult_t,dJw);
  Aut_s =  MortarQuadrature.integrate(f1,Ns_t,Nmult_t,dJw);
  asbMu.localAssembly(umDof,tDof(2:3),-Aut_m);
  asbDu.localAssembly(usDof,tDof(2:3),Aut_s);

  % rhs (jump(eta),t)
  rhsUm(umDof) = rhsUm(umDof) - Aun_m*dTrac(1);
  rhsUs(usDof) = rhsUs(usDof) + Aun_s*dTrac(1);
  rhsUm(umDof) = rhsUm(umDof) - Aut_m*dTrac(2:3);
  rhsUs(usDof) = rhsUs(usDof) + Aut_s*dTrac(2:3);

  % assemble jacobian and rhs of traction balance equations

  % STICK MODE
  if contactState == ContactMode.stick

    asbMt.localAssembly(tDof(1),umDof,-Aun_m');
    asbDt.localAssembly(tDof(1),usDof,Aun_s');
    asbMt.localAssembly(tDof(2:3),umDof,-Aut_m');
    asbDt.localAssembly(tDof(2:3),usDof,Aut_s');

    % (mu_n,g_n)
    rhsT(tDof(1)) = rhsT(tDof(1)) ...
      + MortarQuadrature.integrate(f1,Nmult_n,g_n,dJw);
    % (mu_t,delta_g_t)
    rhsT(tDof(2:3)) = rhsT(tDof(2:3))  ...
      + MortarQuadrature.integrate(f1,Nmult_t,dgt,dJw);
  end

  % SLIP MODE
  if contactState == ContactMode.slip || contactState == ContactMode.newSlip

    % A_nu
    Anu_m = MortarQuadrature.integrate(f1, Nmult_n, Nm_n, dJw);
    Anu_s = MortarQuadrature.integrate(f1, Nmult_n, Ns_n, dJw);
    asbMt.localAssembly(tDof(1),umDof,-Anu_m);
    asbDt.localAssembly(tDof(1),usDof,Anu_s);

    % A_tu (non linear term)
    if obj.state.slip(is) > obj.activeSet.tol.sliding && usDiff > 100*eps
      % compute only on slip terms with sliding large enough
      dtdgt = computeDerTracGap(obj,trac(1),dgt);
      Atu_m = MortarQuadrature.integrate(f2, Nmult_t,dtdgt,Nm_t,dJw);
      Atu_s = MortarQuadrature.integrate(f2, Nmult_t,dtdgt,Ns_t,dJw);
      asbMt.localAssembly(tDof(2:3),umDof,Atu_m);
      asbDt.localAssembly(tDof(2:3),usDof,-Atu_s);

      % A_tn (non linear term)
      dtdtn = computeDerTracTn(obj,is,dgt);
      Atn = MortarQuadrature.integrate(f1,Nmult_t,dtdtn,dJw);
      asbQ.localAssembly(tDof(2:3),tDof(1),-Atn);
    end

    % A_tt
    Att = MortarQuadrature.integrate(f1,Nmult_t,Nmult_t,dJw);
    asbQ.localAssembly(tDof(2:3),tDof(2:3),Att);


    % rhs (mu_n,g_n)
    rhsT(tDof(1)) = rhsT(tDof(1)) + ...
      MortarQuadrature.integrate(f1,Nmult_n,g_n,dJw);

    % rhs -(mu_t,t*_T) (non linear term)
    tT_lim = computeLimitTraction(obj,is,dgt,trac,usDiff);

    %           % check deviation between limit traction and available traction
    %           dev = (R*trac)'*mean(tT_lim,3);
    %           if dev < 0
    %             %trust the limit traction
    %             trac(2:3) = -trac(2:3);
    %             obj.state.traction([3*is-1 3*is]) = -obj.state.traction([3*is-1 3*is]);
    %           end

    % rhs (mu_t,tT) tT and tT_lim in local coordinates (to be consistent with dof definition)
    rhsT(tDof(2:3)) = rhsT(tDof(2:3)) + Att*trac(2:3);
    % reconvert in local coordinates
    %             tT_lim = pagemtimes(R_el,'ctranspose',tT_lim,'none');
    rhsT(tDof(2:3)) = rhsT(tDof(2:3)) - ...
      MortarQuadrature.integrate(f1,Nmult_t,tT_lim,dJw);



    if obj.domains(2).simparams.verbosity > 2
      fprintf('\nelement %i- rhsT: %5.3e %5.3e \n',is,Att*trac(2:3))
      fprintf('element %i- rhsTlim: %5.3e %5.3e \n',is,MortarQuadrature.integrate(f1,Nmult_t,tT_lim,dJw))
      fprintf('------------------------------------ \n')
    end

  end

  % OPEN MODE
  if contactState == ContactMode.open

    % A_oo
    Aoo = MortarQuadrature.integrate(f1,Nmult,Nmult,dJw);
    asbQ.localAssembly(tDof,tDof,Aoo);

    % rhs (mu,t)
    rhsT(tDof) = rhsT(tDof) + ...
      MortarQuadrature.integrate(f1,Nmult,trac,dJw);
  end

end % end inner master elems loop

% assemble matrices into jacobian blocks
obj.addJum(MortarSide.master, asbMu.sparseAssembly());
obj.addJum(MortarSide.slave, asbDu.sparseAssembly());
obj.addJmu(MortarSide.master, asbMt.sparseAssembly());
obj.addJmu(MortarSide.slave, asbDt.sparseAssembly());

obj.Jconstraint = asbQ.sparseAssembly();

obj.addRhs(MortarSide.master,rhsUm);
obj.addRhs(MortarSide.slave,rhsUs);
obj.rhsConstraint = rhsT;

% assemble rhs contribution

% obj.Jmu = asbMu.sparseAssembly();
% obj.Jsu = asbDu.sparseAssembly();
% obj.Jmt = asbMt.sparseAssembly();
% obj.Jst = asbDt.sparseAssembly();
% obj.Jtt = asbQ.sparseAssembly();

end
