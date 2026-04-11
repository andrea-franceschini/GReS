classdef VariablySaturatedFlow < SinglePhaseFlowFVTPFA
  % Variably Saturated flow
  % Subclass of SPFlow
  % Implements Richards equations for unsaturated flow in vadose region
  % It is subclass of SPFlow since most of the methods are the same

  properties
    lw            % upstream mobility
    % JNewt = []      % newton jacobian contribution
    upElem          % upstream elements array for each face
    NLscheme      % newton or picard
  end

  methods (Access = public)
    function obj = VariablySaturatedFlow(domain)
      % initialize as SPFlow class
      obj@SinglePhaseFlowFVTPFA(domain);
    end

    function registerSolver(obj,varargin)

      registerSolver@SinglePhaseFlowFVTPFA(obj,varargin{:});
      % additional logic for richards goes here...
      input = readInput(struct('NLscheme',"newton"),varargin{:});
      obj.NLscheme = input.NLscheme;

    end

    function assembleSystem(obj,dt)

      % compute upstream saturation and mobility
      p = obj.getState(obj.getField());
      [Sw,dSw,d2Sw, obj.lw,dlw] = computeUpElemAndProperties(obj,p);
      % compute obj.H and obj.P
      computeStiffMat(obj,obj.lw);
      computeCapMat(obj,Sw,dSw);

      ents = obj.domain.dofm.getActiveEntities(obj.fieldId);


      if obj.steadyState
        rhs = obj.H*p(ents);
        J = obj.H;
      else
        pOld = obj.getStateOld(obj.getField());
        rhs = obj.H*p(ents) + (obj.P/dt)*(p(ents) - pOld(ents));
        J = obj.H + obj.P/dt;
      end


      if isNewtonNLSolver(obj)
        Jh = computeJacobianPartJhNewton(obj,p,dlw);
        Jp = computeJacobianPartJpNewton(obj,dt,p,pOld,Sw,dSw,d2Sw);
        J = J + Jh;
        if ~obj.steadyState
          J = J + Jp/dt;
        end
      end

      obj.domain.J{obj.fieldId,obj.fieldId} = J;

      gamma = obj.domain.materials.getFluid().getSpecificWeight();
      if gamma > 0
        % add rhs gravity contribution
        obj.domain.rhs{obj.fieldId} = rhs + getRhsGravity(obj,obj.lw);
      else
        obj.domain.rhs{obj.fieldId} = rhs;
      end


    end



    function states = finalizeState(obj,p,t)
      % Compute the posprocessing variables for the module.
      gamma = obj.domain.materials.getFluid().getSpecificWeight();
      if gamma>0
        zbc = obj.grid.cells.center(:,3);
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



    function gTerm = getRhsGravity(obj,lw)

      dofm = obj.domain.dofm;
      nCells = dofm.getNumbDoF(obj.fieldId);
      neigh = obj.grid.faces.neighbors(obj.isIntFaces,:);
      gTerm = accumarray(neigh(:),[lw.*obj.rhsGrav; ...
        -lw.*obj.rhsGrav],[nCells,1]);

      gTerm = gTerm(dofm.getActiveEntities(obj.fieldId));
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

    function writeSolution(obj,fac,tID)
      satOld = getStateOld(obj,"saturation");
      satCurr = getState(obj,"saturation");
      pOld = getStateOld(obj,"pressure");
      pCurr = getState(obj,"pressure");

      obj.domain.outstate.results(tID).pressure = pCurr*fac+pOld*(1-fac);
      obj.domain.outstate.results(tID).saturation = satCurr*fac+satOld*(1-fac);
    end

    function [ents,vals] = getBC(obj,bcId,t)
      % getBC - function to find the value and the location for the
      % boundary condition.
      
      % overrides the SinglePhaseFlowFVTPFA getBC() to add the newton
      % contribution to the boundary condition values

      faces = obj.grid.faces;
      cells = obj.grid.cells;

      bc = obj.domain.bcs;
      mat = obj.domain.materials;

      bcFld = bc.getField(bcId);
      type = bc.getType(bcId);

      if strcmp(type,'seepage')
        assert(isEssential(bc,bcId),"Boundary condition of type 'seepage' must be essential")
      end

      if bcFld == entityField.surface && isEssential(bc,bcId)

        faceId = bc.getSourceEntities(bcId);
        zf = faces.center(faceId,3);
        srcVal = bc.getSourceVals(bcId,t);
        gamma = mat.getFluid().getSpecificWeight();

        if strcmp(type,"seepage")
          srcVal = gamma*(srcVal-zf);
          srcVal(srcVal<=0)=0.;
        end

        ents = sum(faces.neighbors(faceId,:),2);
        p = getState(obj,obj.getField());

        %
        [mob, dmob] = obj.computeMobilityBoundary( ...
          obj.domain.state.data.pressure(ents),srcVal,faceId);
        tr = obj.trans(faceId);
        dz = cells.center(ents,3) - zf;
        dirJ = mob.*tr;
        potential = (p(ents) - srcVal) + gamma*dz;
        q = dirJ.*potential;
        if isNewtonNLSolver(obj)
          % Contibution of Jh part in the boundary
          dirJh = dmob.*tr;
          dirJ = dirJ + dirJh.*potential;
        end
        vals = [dirJ,q]; % {JacobianVal,rhsVal]

      else

        ents = bc.getTargetEntities(bcId);
        vals = bc.getVals(bcId,t);

      end

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
    function [Sw,dSw,d2Sw] = computeSaturation(obj,p)
      % COMPUTESATURATION compute the saturation and it's derivatives
      cells = obj.grid.cells;
      Sw = zeros(cells.num,1);
      dSw = zeros(cells.num,1);
      d2Sw = zeros(cells.num,1);

      for m = 1:cells.nTag
        isElMat = cells.tag == m;
        p = p(isElMat);
        mat = obj.domain.materials.getMaterial(m);
        Sws = mat.PorousRock.getMaxSaturation();
        Swr = mat.PorousRock.getResidualSaturation();
        [Sw(isElMat), dSw(isElMat), d2Sw(isElMat)] = mat.PorousRock.Curves.computeSwAnddSw(p);
        Sw(isElMat) = Swr + (Sws-Swr)*Sw(isElMat);
        dSw(isElMat) = (Sws-Swr)*dSw(isElMat);
        d2Sw(isElMat) = (Sws-Swr)*d2Sw(isElMat);
      end
    end

    function [lw,dlw] = computeMobility(obj,p)
      % COMPUTEMOBILITY compute the mobility and its derivatives
      % for the upstream elements for each face
      nIntFaces = length(obj.upElem);
      lw = zeros(nIntFaces,1);
      dlw = zeros(nIntFaces,1);
      cells = obj.grid.cells;
      matUpElem = cells.tag(obj.upElem);
      mat = obj.domain.materials;
      for m = 1:cells.nTag
        isElMat = matUpElem == m;
        p = p(obj.upElem(isElMat));
        [lw(isElMat), dlw(isElMat)] = mat.getMaterial(m).PorousRock.Curves.computeRelativePermeability(p);
        % [lw(isElMat), dlw(isElMat)] = mat.getMaterial(m).RelativePermCurve.interpTable(p);
      end
      mu = mat.getFluid().getDynViscosity();
      lw = lw/mu;
      dlw = dlw/mu;
    end

    function [lpt, dlpt] = computeMobilityBoundary(obj,pcells,pface,faceId)
      % COMPUTEMOBILITYBOUNDARY compute the mobility for the
      % upstream elements in the boundary
      faces = obj.grid.faces;
      cells = obj.grid.cells;

      elms = faces.neighbors(faceId,:);
      elms = elms(elms~=0);
      materialsID = cells.tag(elms);
      mat = obj.domain.materials;

      % Find the direction of the flux;
      gamma = mat.getFluid().getSpecificWeight();
      if gamma > 0
        zfaces = faces.center(faceId,3);
        cellz = cells.center(elms,3);
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
        [krel,dkrel] = mat.getMaterial(materialsID(i)).PorousRock.Curves.computeRelativePermeability(pres(i));
        lpt(i) = krel/mu;
        dlpt(i) = -dkrel/mu;
      end
      dlpt = dlpt.*lElemIsUp;
    end

    function [Sw,dSw,d2Sw,lw,dlw] = computeUpElemAndProperties(obj,p)
      % compute upstream elements for each face
      % interpolate effective saturation and relative permeability
      % and first derivatives
      % compute also second derivative for saturation
      neigh = obj.grid.faces.neighbors(obj.isIntFaces,:);
      gamma = obj.domain.materials.getFluid().getSpecificWeight();
      if gamma > 0
        zVec = obj.grid.cells.center(:,3);
        zNeigh = zVec(neigh);
        lElemIsUp = (p(neigh(:,1)) - p(neigh(:,2))) + gamma*(zNeigh(:,1) - zNeigh(:,2)) >= 0;
      else
        lElemIsUp = p(neigh(:,1)) >= p(neigh(:,2));
      end
      obj.upElem(lElemIsUp) = neigh(lElemIsUp,1);
      obj.upElem(~lElemIsUp) = neigh(~lElemIsUp,2);
      [Sw,dSw,d2Sw] = computeSaturation(obj,p);
      dSw = - dSw;
      [lw,dlw] = computeMobility(obj,p);
      dlw = - dlw;
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

      gamma = obj.domain.materials.getFluid.getSpecificWeight();
      % subCells = obj.domain.dofm.getFieldCells(obj.getField());

      % Get pairs of faces that contribute to the subdomain
      neigh = obj.grid.faces.neighbors(obj.isIntFaces,:);
      zVec = obj.grid.cells.center(:,3);
      zNeigh = zVec(neigh);
      pTau = pTau(neigh);
      nneigh = length(dMob);
      % [~,~,reorder] = unique([neigh(:,1); neigh(:,2); obj.upElem'; subCells]);
      [~,~,reorder] = unique([neigh(:,1); neigh(:,2); obj.upElem']);
      neigh1 = reorder(1:nneigh);
      neigh2 = reorder(nneigh+1:2*nneigh);
      neighU = reorder(2*nneigh+1:3*nneigh);

      % Transmissibility of internal faces
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
      mat = obj.domain.materials;
      beta = mat.getFluid().getFluidCompressibility();
      subCells = obj.domain.dofm.getFieldCells(obj.getField());
      nSubCells = length(subCells);
      poroMat = zeros(nSubCells,1);
      alphaMat = zeros(nSubCells,1);
      cells = obj.grid.cells;
      for m = 1:cells.nTag
        alphaMat(m) = getRockCompressibility(obj,m);
        poroMat(m) = mat.getMaterial(m).PorousRock.getPorosity();
      end

      alphaMat = alphaMat(cells.tag(subCells));
      poroMat = poroMat(cells.tag(subCells));
 
      % Define some parameters.
      pTmp = pTmp(subCells);
      pOld = pOld(subCells);
      pdiff = pTmp-pOld;

      % Computing the Jacobian Part.
      Jp = beta*alphaMat.*Stau + (2*alphaMat+beta*poroMat).*dStau + poroMat.*ddStau;
      Jp = cells.volume(subCells).*Jp.*pdiff;

      nDoF = obj.domain.dofm.getNumbDoF(obj.getField());
      [~,~,dof] = unique(subCells);
      Jp = sparse(dof,dof,Jp,nDoF,nDoF);
    end


    function out = isNewtonNLSolver(obj)
      out = strcmp("newton",obj.NLscheme);
    end

  end

  methods (Static)
    

    % function out = getField()
    %   out = "pressure";
    % end

  end
end

