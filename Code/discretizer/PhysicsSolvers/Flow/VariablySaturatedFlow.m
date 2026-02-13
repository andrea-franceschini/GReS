classdef VariablySaturatedFlow < SinglePhaseFlowFVTPFA
  % Variably Saturated flow
  % Subclass of SPFlow
  % Implements Richards equations for unsaturated flow in vadose region
  % It is subclass of SPFlow since most of the methods are the same

  properties
    lwkpt           % mobility
    % JNewt = []      % newton jacobian contribution
    upElem          % upstream elements array for each face
  end

  methods (Access = public)
    function obj = VariablySaturatedFlow(domain)
      % initialize as SPFlow class
      obj@SinglePhaseFlowFVTPFA(domain);
    end

    function registerSolver(obj,solverInput)
      registerSolver@SinglePhaseFlowFVTPFA(obj,solverInput);
      % additional logic for richards goes here...
    end

    function states = finalizeState(obj,p,t)
      % Compute the posprocessing variables for the module.
      gamma = obj.domain.materials.getFluid().getSpecificWeight();
      if gamma>0
        zbc = obj.mesh.cellCentroid(:,3);
        states.potential = p + gamma*zbc;
        states.head = zbc+p/gamma;
      end
      [mob ,~] = computeMobility(obj,p);
      states.flux = computeFlux(obj,p,mob,t);
      states.perm = printPermeab(obj);
      states.saturation = computeSaturation(obj,p);
      states.pressure = p;
      % states.mass = checkMassCons(obj,mob,potential);
    end

    function J = computeMat(obj,dt)
      theta = obj.domain.simparams.theta;

      % get pressure state
      p = getState(obj,obj.getField());
      pOld = getStateOld(obj,obj.getField());
      pkpt = theta*p + (1 - theta)*pOld;
      
      [Swkpt,dSwkpt,d2Swkpt, obj.lwkpt,dlwkpt] = computeUpElemAndProperties(obj,pkpt);
      computeStiffMat(obj,obj.lwkpt);
      computeCapMat(obj,Swkpt,dSwkpt);
      J = theta*obj.H + obj.P/dt;
      % obj.H = theta*obj.H;
      % obj.P = theta*obj.P/dt;
      if isNewtonNLSolver(obj.domain.simparams)
        J = J + computeJacobianPartJhNewton(obj,pkpt,dlwkpt);
        J = J + computeJacobianPartJpNewton(obj,dt,p,pOld,Swkpt,dSwkpt,d2Swkpt);

        % obj.H = obj.H + computeJacobianPartJhNewton(obj,pkpt,dlwkpt);
        % obj.P = obj.P + computeJacobianPartJpNewton(obj,dt,obj.domain.state.data.pressure, ...
        %   obj.domain.stateOld.data.pressure,Swkpt,dSwkpt,d2Swkpt);
      end
      % J = obj.H + obj.P;
    end

    function rhs = computeRhs(obj,dt)
      % Compute the residual of the flow problem

      % get pressure state
      p = getState(obj,obj.getField());
      pOld = getStateOld(obj,obj.getField());

      theta = obj.domain.simparams.theta;
      ents = obj.domain.dofm.getActiveEntities(obj.getField());

      % computing the rhs
      rhs = obj.H*(theta*p(ents) + (1-theta)*pOld(ents));
      rhs = rhs + (obj.P/dt)*(p(ents) - pOld(ents));

      % adding gravity contribution
      gamma = obj.domain.materials.getFluid().getSpecificWeight();
      if gamma > 0
        rhs = rhs + finalizeRHSGravTerm(obj,obj.lwkpt);
      end
    end

    function updateState(obj,dSol)
      dofm = obj.domain.dofm;
      if nargin > 1
        ents = dofm.getActiveEntities(obj.getField());
        state = getState(obj);
        p = state.data.pressure;
        state.data.pressure(ents) = p(ents) + dSol(dofm.getDoF(obj.fieldId));
        state.data.saturation = computeSaturation(obj,p(ents));
      end
    end

    function writeMatFile(obj,fac,tID)
      satOld = getStateOld(obj,"saturation");
      satCurr = getState(obj,"saturation");
      pOld = getStateOld(obj,"pressure");
      pCurr = getState(obj,"pressure");

      obj.domain.outstate.matFile(tID).pressure = pCurr*fac+pOld*(1-fac);
      obj.domain.outstate.matFile(tID).saturation = satCurr*fac+satOld*(1-fac);
    end

    function [ents,vals] = getBC(obj,id,t)

      % update to new logic. can we reuse entirely the SinglePhaseFlow
      
      mat = obj.domain.materials;
      bc = obj.domain.bcs;

      switch bc.getCond(id)
        case {'NodeBC','ElementBC'}
          ents = bc.getEntities(id);
          vals = bc.getVals(id,t);
        case 'SurfBC'
          % faceID = bc.getEntities(id);
          % ents = sum(obj.faces.faceNeighbors(faceID,:),2);
          % v = bc.getVals(id,t);
          % [ents,~,ind] = unique(ents);

          [faceID, faceOrder] = sort(bc.getEntities(id));
          ents = sum(obj.faces.faceNeighbors(faceID,:),2);
          v(faceOrder,1) = bc.getVals(id,t);
          switch bc.getType(id)
            case 'Neumann'
              vals = vecnorm(obj.faces.faceNormal(faceID,:),2,2).*v;
              % area = vecnorm(obj.faces.faceNormal(faceID,:),2,2).*v;
              % vals = accumarray(ind, area);
            case 'Dirichlet'
              % theta = obj.simParams.theta;
              gamma = mat.getFluid().getSpecificWeight();
              [mob, dmob] = obj.computeMobilityBoundary( ...
                obj.domain.state.data.pressure(ents),v,faceID);
              tr = obj.trans(faceID);
              press = obj.domain.state.data.pressure(ents) - v;
              gravT =  gamma*(obj.mesh.cellCentroid(ents,3) ...
                - obj.faces.faceCentroid(faceID,3));
              dirJ = mob.*tr;
              q = dirJ.*(press+gravT);
              if isNewtonNLSolver(obj.domain.simparams)
                % Contibution of Jh part in the boundary
                dirJh = dmob.*tr;
                dirJ = dirJ + dirJh.*(press+gravT);
              end
              vals = [dirJ,q]; % {JacobianVal,rhsVal]
              % vals = [dirJ,accumarray(ind,q)]; % {JacobianVal,rhsVal]
            case 'Seepage'
              gamma = mat.getFluid().getSpecificWeight();
              assert(gamma>0.,'To impose Seepage boundary condition is necessary the fluid specify weight be bigger than zero!');

              % theta = obj.simParams.theta;
              zbc = obj.faces.faceCentroid(faceID,3);
              href = v;
              v = gamma*(href(1)-zbc);

              v(v<=0)=0.;
              [mob, dmob] = obj.computeMobilityBoundary( ...
                obj.domain.state.data.pressure(ents),v,faceID);
              tr = obj.trans(faceID);
              press = obj.domain.state.data.pressure(ents) - v;
              gravT = gamma*(obj.mesh.cellCentroid(ents,3) ...
                - obj.faces.faceCentroid(faceID,3));
              dirJ = mob.*tr;
              q = dirJ.*(press+gravT);
              if isNewtonNLSolver(obj.domain.simparams)
                % Contibution of Jh part in the boundary
                dirJh = dmob.*tr;
                dirJ = dirJ + dirJh.*(press+gravT);
              end
              vals = [dirJ,q]; % {JacobianVal,rhsVal]

              % % pos=v>=0;
              % % % v(v<=0)=0;  % Atmosferic pressure
              % % % resize the number of boundary condition.
              % % [ents,~,ind] = unique(ents(pos));
              % % v=v(pos); faceID = faceID(pos);
              % % [mob, dmob] = obj.computeMobilityBoundary(state.pressure(ents),v,faceID);
              % % tr = obj.trans(faceID);
              % % press = state.pressure(ents) - v;
              % % gravT = gamma*(obj.elements.cellCentroid(ents,3) ...
              % %    - obj.faces.faceCentroid(faceID,3));
              % % dirJ = mob.*tr;
              % % q = dirJ.*(press+gravT);
              % % if isNewtonNLSolver(obj.domain.simparams)
              % %    % Contibution of Jh part in the boundary
              % %    dirJh = dmob.*tr.*press;
              % %    dirJ = dirJ + dirJh;
              % % end
              % % vals = [dirJ,accumarray(ind,q)]; % {JacobianVal,rhsVal]
              % % % vals = [dirJ,q]; % {JacobianVal,rhsVal]
          end

        case 'VolumeForce'
          v = bc.getVals(id,t);
          ents = bc.getEntities(id);
          vals = v.*obj.mesh.cellVolume(ents);
      end
      % get local dof numbering
      %dof = obj.domain.dofm.getLocalDoF(obj.fieldId,ents);
    end

    function [cellStr,pointStr] = buildPrintStruct(obj,state)
      pointStr = repmat(struct('name', 1, 'data', 1), 1, 1);
      pointStr(1).name = 'flux';
      pointStr(1).data = state.flux;

      cellStr = repmat(struct('name', 1, 'data', 1), 3, 1);
      cellStr(1).name = 'pressure';
      cellStr(1).data = state.pressure;      
      cellStr(2).name = 'saturation';
      cellStr(2).data = state.saturation;
      cellStr(3).name = 'permeability';
      cellStr(3).data = state.perm;
      if isfield(state,"potential")
        cellStr(4).name = 'potential';
        cellStr(4).data = state.potential;      
        cellStr(5).name = 'piezometric head';
        cellStr(5).data = state.head;
      end
      % cellStr(6).name = 'mass_cons';
      % cellStr(6).data = state.mass;
    end

    function out = isFEM(obj)
      out = false;
    end

    function out = isTPFA(obj)
      out = true;
    end

    function str = typeDiscretization(obj)
      str = "FVTPFA";
    end

    function out = isLinear(obj)
      out = false;
    end

  end

  methods (Access = private)
    function [Swkpt,dSwkpt,d2Swkpt] = computeSaturation(obj,pkpt)
      % COMPUTESATURATION compute the saturation and it's derivatives
      Swkpt = zeros(obj.mesh.nCells,1);
      dSwkpt = zeros(obj.mesh.nCells,1);
      d2Swkpt = zeros(obj.mesh.nCells,1);
      for m = 1:obj.mesh.nCellTag
        isElMat = obj.mesh.cellTag == m;
        p = pkpt(isElMat);
        mat = obj.domain.materials.getMaterial(m);
        Sws = mat.PorousRock.getMaxSaturation();
        Swr = mat.PorousRock.getResidualSaturation();
        [Swkpt(isElMat), dSwkpt(isElMat), d2Swkpt(isElMat)] = mat.Curves.computeSwAnddSw(p);
        Swkpt(isElMat) = Swr + (Sws-Swr)*Swkpt(isElMat);
        dSwkpt(isElMat) = (Sws-Swr)*dSwkpt(isElMat);
        d2Swkpt(isElMat) = (Sws-Swr)*d2Swkpt(isElMat);
      end
    end

    function [lwkpt,dlwkpt] = computeMobility(obj,pkpt)
      % COMPUTEMOBILITY compute the mobility and it's derivatives
      % for the upstream elements for each face
      nIntFaces = length(obj.upElem);
      lwkpt = zeros(nIntFaces,1);
      dlwkpt = zeros(nIntFaces,1);
      matUpElem = obj.mesh.cellTag(obj.upElem);
      mat = obj.domain.materials;
      for m = 1:obj.mesh.nCellTag
        isElMat = matUpElem == m;
        p = pkpt(obj.upElem(isElMat));
        [lwkpt(isElMat), dlwkpt(isElMat)] = mat.getMaterial(m).Curves.computeRelativePermeability(p);
        % [lwkpt(isElMat), dlwkpt(isElMat)] = mat.getMaterial(m).RelativePermCurve.interpTable(p);
      end
      mu = mat.getFluid().getDynViscosity();
      lwkpt = lwkpt/mu;
      dlwkpt = dlwkpt/mu;
    end

    function [lpt, dlpt] = computeMobilityBoundary(obj,pcells,pface,faceID)
      % COMPUTEMOBILITYBOUNDARY compute the mobility for the
      % upstream elements in the boundary
      elms = obj.faces.faceNeighbors(faceID,:);
      elms = elms(elms~=0);
      materialsID = obj.mesh.cellTag(elms);
      mat = obj.domain.materials;

      % Find the direction of the flux;
      gamma = mat.getFluid().getSpecificWeight();
      if gamma > 0
        zfaces = obj.faces.faceCentroid(faceID,3);
        cellz = obj.mesh.cellCentroid(elms,3);
        lElemIsUp = (pcells - pface) + gamma*(cellz- zfaces) >= 0;
      else
        lElemIsUp = pcells >= pface;
      end

      % Find the upstream pressure.
      pres = zeros(length(lElemIsUp),1);
      pres(lElemIsUp) = pcells(lElemIsUp);
      pres(~lElemIsUp) = pface(~lElemIsUp);

      mu = mat.getFluid().getDynViscosity();
      lpt = zeros(length(pcells),1);
      dlpt = zeros(length(pcells),1);
      for i=1:length(materialsID)
        % [krel,dkrel] = mat.computeRelativePermeability(pres(i),materials(i));
        [krel,dkrel] = mat.getMaterial(materialsID(i)).Curves.computeRelativePermeability(pres(i));
        lpt(i) = krel/mu;
        dlpt(i) = -dkrel/mu;
      end
      dlpt = dlpt.*lElemIsUp;
    end

    function [Swkpt,dSwkpt,d2Swkpt,lwkpt,dlwkpt] = computeUpElemAndProperties(obj,pkpt)
      % compute upstream elements for each face
      % interpolate effective saturation and relative permeability
      % and first derivatives
      % compute also second derivative for saturation
      neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
      gamma = obj.domain.materials.getFluid().getSpecificWeight();
      if gamma > 0
        zVec = obj.elements.mesh.cellCentroid(:,3);
        zNeigh = zVec(neigh);
        lElemIsUp = (pkpt(neigh(:,1)) - pkpt(neigh(:,2))) + gamma*(zNeigh(:,1) - zNeigh(:,2)) >= 0;
      else
        lElemIsUp = pkpt(neigh(:,1)) >= pkpt(neigh(:,2));
      end
      obj.upElem(lElemIsUp) = neigh(lElemIsUp,1);
      obj.upElem(~lElemIsUp) = neigh(~lElemIsUp,2);
      [Swkpt,dSwkpt,d2Swkpt] = computeSaturation(obj,pkpt);
      dSwkpt = - dSwkpt;
      [lwkpt,dlwkpt] = computeMobility(obj,pkpt);
      dlwkpt = - dlwkpt;
    end

    function Jh = computeJacobianPartJhNewton(obj,pTau,dMob)
      % COMPUTEJACOBIANFORNEWTONJFPART Method to compute the Jh part
      % of the jacobian used in the Newton-Raphson.
      %% Equation for the part.
      % $$ J_{h} = \theta p_{\tau}^{m}
      % \frac{\partial\lambda_{u}^{t+dt,m}}{\partial p_{t+dt}} \mathbf{T}$$
      % $$ J_{f} = \theta \gamma
      % \frac{\partial\lambda_{u}^{t+dt,m}}{\partial p_{t+dt}} \mathbf{T}$$
      % $$ J_{hf} = J_{h} + J_{f}$$

      % Define some values.
      % theta = obj.simParams.theta;
      gamma = obj.domain.materials.getFluid.getSpecificWeight();
      % subCells = obj.domain.dofm.getFieldCells(obj.getField());

      % Get pairs of faces that contribute to the subdomain
      neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
      zVec = obj.mesh.cellCentroid(:,3);
      zNeigh = zVec(neigh);
      pTau = pTau(neigh);
      nneigh = length(dMob);
      % [~,~,reorder] = unique([neigh(:,1); neigh(:,2); obj.upElem'; subCells]);
      [~,~,reorder] = unique([neigh(:,1); neigh(:,2); obj.upElem']);
      neigh1 = reorder(1:nneigh);
      neigh2 = reorder(nneigh+1:2*nneigh);
      neighU = reorder(2*nneigh+1:3*nneigh);

      % Transmissibility of internal faces
      % Jh = theta*dMob.*obj.trans(obj.isIntFaces);
      Jh = dMob.*obj.trans(obj.isIntFaces);
      Jh = Jh.*(pTau(:,1) - pTau(:,2) + gamma*(zNeigh(:,1) - zNeigh(:,2)));

      % Computing the Jacobian Part.
      nDoF = obj.domain.dofm.getNumbDoF(obj.getField());
      Jh = sparse( [neigh1; neigh2], [neighU; neighU], ...
        [Jh; -Jh], nDoF, nDoF);
    end

    function Jp = computeJacobianPartJpNewton(obj,dt,pTmp,pOld,Stau,dStau,ddStau)
      % COMPUTEJACOBIANFORNEWTONPPART Method to compute the Jp part
      % of the jacobian used in the Newton-Raphson.
      %% Equation for the part.
      % $$ J_{p} = \left(\theta/\delta t\right) \left[
      % \alpha^{e}\beta S_{\tau}^{e,m}
      % + \left(2\alpha^{e}+\beta\phi_{\tau}^{e}\right)\frac{\partial S_{\tau}^{e,m}}{\partial p}
      % + \phi_{\tau}^{e} \frac{\partial^{2} S_{t+dt}^{m,e}}{\partial p_{t+dt}^{2}} \right]
      % \left(\frac{p_{t+dt}^{m}-p_{t}}{dt}\right)
      % \left|\Omega\right|^{e}$$

      % Define some constant.
      % theta = obj.simParams.theta;
      mat = obj.domain.materials;
      beta = mat.getFluid().getFluidCompressibility();
      subCells = obj.domain.dofm.getFieldCells(obj.getField());
      nSubCells = length(subCells);
      poroMat = zeros(nSubCells,1);
      alphaMat = zeros(nSubCells,1);
      for m = 1:obj.mesh.nCellTag
        alphaMat(m) = getRockCompressibility(obj,obj.mesh.cellTag(m));
        % if ~ismember(m,obj.domain.dofm.getFieldCellTags({obj.getField(),"displacements"}))
        %   % compute alpha only if there's no coupling in the subdomain
        %   alphaMat(m) = mat.getMaterial(m).ConstLaw.getRockCompressibility();
        % end
        poroMat(m) = mat.getMaterial(m).PorousRock.getPorosity();
      end

      alphaMat=alphaMat(obj.mesh.cellTag(subCells));
      poroMat=poroMat(obj.mesh.cellTag(subCells));

      % Define some parameters.
      pTmp = pTmp(subCells);
      pOld = pOld(subCells);
      pdiff = pTmp-pOld;

      % Computing the Jacobian Part.
      % Jp = poroMat.*ddStau;
      % Jp = (theta/dt)*obj.elements.vol(subCells).*Jp.*pdiff;
      Jp = beta*alphaMat.*Stau + (2*alphaMat+beta*poroMat).*dStau + poroMat.*ddStau;
      Jp = (1./dt)*obj.elements.mesh.cellVolume(subCells).*Jp.*pdiff;

      nDoF = obj.domain.dofm.getNumbDoF(obj.getField());
      [~,~,dof] = unique(subCells);
      Jp = sparse(dof,dof,Jp,nDoF,nDoF);
    end

  end

  methods (Static)
    

    % function out = getField()
    %   out = "pressure";
    % end

  end
end

