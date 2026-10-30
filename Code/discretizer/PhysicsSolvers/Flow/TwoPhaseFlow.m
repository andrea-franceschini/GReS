classdef TwoPhaseFlow < SinglePhaseFlowFVTPFA
  % Immiscible Two-phase flow solver using MRST ADI framework
  % we assume that phases are immiscible and capillary effects are
  % negligible

  % Uses TPFA operators

  properties

    ops                                               % struct collecting AD operators
    props = struct('wet',[],'nwet',[],'rock',[])      % struct collecting fluid props

  end

  methods
    function obj = TwoPhaseFlow(domain)

      obj@SinglePhaseFlowFVTPFA(domain);

    end


    function registerSolver(obj,varargin)


      registerSolver@SinglePhaseFlowFVTPFA(obj,varargin{:});
      dofm = obj.domain.dofm;

      regions = getTargetRegions(dofm,"pressure");
      dofm.registerVariable("saturation",entityField.cell,1,regions);
      state = getState(obj);
      state.saturation = zeros(numel(state.pressure),1);
      setState(obj,state);

      % set ADI operators
      G = obj.grid;
      Nglobal = G.faces.neighbors(obj.isIntFaces,:);
      cellIdx = obj.domain.dofm.getFieldCells("pressure");
      nc = numel(cellIdx);
      nf = size(Nglobal,1);
      [found,N] = ismember(Nglobal,cellIdx);
      assert(all(found(:)), "Internal faces must connect active pressure cells.");

      D = sparse([(1:nf)'; (1:nf)'], [N(:,1); N(:,2)], ...
        [-ones(nf,1); ones(nf,1)], nf, nc);
      obj.ops.grad  = @(x)  D*x;
      obj.ops.div   = @(x) -D'*x;
      obj.ops.avg   = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));
      obj.ops.upw   = @(flag, x) flag.*x(N(:, 1)) + ~flag.*x(N(:, 2));
      obj.ops.gradz = obj.ops.grad(G.cells.center(cellIdx,3));


    end


    function setFluidProperties(obj,input)

      fluidParams = struct('viscosity',[],'compressibility',0.0,'relPerm',[],'density',[]);
      rockParams = struct('refPressure',missing);
      default = struct('wet',fluidParams,'nwet',fluidParams,'rock',rockParams);

      input = readInput(default,input);

      obj.props.wet = input.wet;
      obj.props.nwet = input.nwet;

      % rock properties
      dofm = obj.domain.dofm;
      regions = getTargetRegions(dofm,"pressure");
      nc = dofm.getNumbDoF("pressure");
      pv0 = zeros(nc,1);
      cr = zeros(nc,1);
      g = obj.domain.grid;
      mat = obj.domain.materials;


      for i = regions'
        rock = mat.getPorousRock(i);
        cIdx = g.cells.tag == i;
        poro0 = rock.getPorosity;
        pv0(cIdx) = g.cells.volume(cIdx).*poro0;
        cr(cIdx) = rock.getCompressibility;
      end

      % pore volume function
      pref = input.rock.refPressure;
      if ~ismissing(pref)
        pV   = @(p) pv0 .* exp(cr .* (p - pref) );
      else
        pV = @(p) 0*p + pv0;
      end

      obj.props.rock.poreVolume = pV;

    end


    function applyBC(obj,bcId,t)

      if ~BCapplies(obj,bcId)
        return
      end

      bcs = obj.domain.boundaries;
      type = bcs.getType(bcId);

      if contains(type,"well")
        % well equations


      else
        applyBC@SinglePhaseFlowFVTPFA(obj,bcId,t);
      end



    end


    function applyDirVal(obj,bcId,varargin)
      if ~BCapplies(obj,bcId)
        return
      end

      if ~isEssential(obj.domain.bcs,bcId)
        return
      end

      if nargin == 3
        type = obj.domain.bcs.getType(bcId);
        assert(~BCtype.isCustomBC(type),"Base applyBC method for custom BC requires " + ...
          "dofs and values to be already specified")
        t = varargin{1};
        dofs = obj.domain.bcs.getStateDofs(bcId);
        vals = obj.domain.bcs.getVals(bcId,t);
      elseif nargin==4
        dofs = varargin{1};
        vals = varargin{2};
      else
        error("Wrong number of input for applyBC function")
      end

      bcVar = obj.domain.bcs.getVariable(bcId);

      v = obj.getState(bcVar);
      v(dofs) = vals;
      obj.setState(v,bcVar);

      % apply Dirichlet condition also to ADI variable!
      if strcmp(bcVar,"pressure")
        obj.ops.p.val(dofs) = vals;
      elseif strcmp(bcVar,"saturation")
        obj.ops.sw.val(dofs) = vals;
      end

    end


    function initialize(obj)

      % material input is provisionally here
      % units SI

      if isempty(obj.props.wet)

        error('Missing fluid properties. Fluid properties must be passed using the "setFluidParameters()" method');

      end

      % vertical equilibrium
      % g = 9.8066;     % gravity acceleration
      % equil = ode23(@(z,p) g.* rhoO(p), [0, max(G.cells.centroids(:,3))], pR);
      cIdx = obj.domain.dofm.getActiveEntities('pressure');
      p0    = obj.getState.pressure(cIdx);
      sW0   = obj.getState.saturation(cIdx);
      [obj.ops.p, obj.ops.sw]    = initVariablesADI(p0, sW0);

      gravity reset off
      msg = sprintf("Gravity not yet supported in %s solver",class(obj));
      gresLog().log(3,msg);


    end



    function assembleSystem(obj,dt)

      cells = obj.domain.dofm.getFieldCells("pressure");
      curr  = getState(obj);
      prev = getStateOld(obj);

      pOld  = prev.pressure(cells);
      sWOld = prev.saturation(cells);

      [p,sW] = initVariablesADI( ...
        curr.pressure(cells),curr.saturation(cells));

      pv = obj.props.rock.poreVolume;

      [wet,nwet] = obj.computeFlowResidual( ...
        dt,p,sW,pOld,sWOld,pv(p),pv(pOld));


      % residual equations
      dofm = obj.domain.dofm;
      eqIdx = [dofm.getVariableId("pressure"); ...
        dofm.getVariableId("saturation")];

      obj.domain.rhs{eqIdx(1)} = wet.val;
      obj.domain.rhs{eqIdx(2)} = nwet.val;
      obj.domain.J(eqIdx(1),eqIdx) = wet.jac;
      obj.domain.J(eqIdx(2),eqIdx) = nwet.jac;

    end


    function [cellData,pointData] = writeVTK(obj,fac,t)

      p = obj.domain.state.interpolate(fac,"pressure");
      sat = obj.domain.state.interpolate(fac,"saturation");

      outPrint = finalizeState(obj,p,t);

      [cellData,pointData] = buildPrintStruct(obj,outPrint);

      cellData(end+1).name = 'saturation';
      cellData(end).data = sat;

    end

    function writeSolution(obj,fac,tID)
      p = obj.domain.state.interpolate(fac,"pressure");
      sat = obj.domain.state.interpolate(fac,"saturation");
      obj.domain.outstate.results(tID).pressure = p;
      obj.domain.outstate.results(tID).saturation = sat;
    end



    function updateState(obj,du)


      % update the adi variable as well as the state object
      dofm = obj.domain.dofm;
      dp = du(dofm.getDoF(dofm.getVariableId("pressure")));
      ds = du(dofm.getDoF(dofm.getVariableId("saturation")));
      [p,sw] = deal(obj.ops.p,obj.ops.sw);
      assert(all(isfinite(dp(:))) && all(isfinite(ds(:))), ...
        "Nonfinite Newton correction: stop before corrupting the state.");
      p.val = p.val + dp;
      sw.val = sw.val + ds;
      % regularize saturation
      sw.val = max( min(sw.val, 1), 0);
      [obj.ops.p,obj.ops.sw] = deal(p,sw);

      setState(obj,value(p),"pressure");
      setState(obj,value(sw),"saturation");

    end



    function [wet,nwet] = computeFlowResidual(obj,dt,p,sW,pOld,sWOld,vol,volOld)

      % Preserve these for helpers that currently access obj.ops.p.
      obj.ops.p  = p;
      obj.ops.sw = sW;

      T = obj.trans(obj.isIntFaces);

      w  = obj.props.wet;
      nw = obj.props.nwet;

      rW   = w.density(p);
      rWOld  = w.density(pOld);
      rNW  = nw.density(p);
      rNWOld = nw.density(pOld);

      mobW  = w.relPerm(sW)./w.viscosity;
      mobNW = nw.relPerm(1-sW)./nw.viscosity;

      g = norm(gravity);
      dp = obj.ops.grad(p);

      dpW  = dp + g*obj.ops.avg(rW).*obj.ops.gradz;
      dpNW = dp + g*obj.ops.avg(rNW).*obj.ops.gradz;

      vW = -obj.ops.upw(value(dpW)<=0,rW.*mobW).*T.*dpW;
      vNW = -obj.ops.upw(value(dpNW)<=0,rNW.*mobNW).*T.*dpNW;

      wet = (vol.*rW.*sW - volOld.*rWOld.*sWOld)/dt ...
        + obj.ops.div(vW);

      nwet = (vol.*rNW.*(1-sW) - volOld.*rNWOld.*(1-sWOld))/dt ...
        + obj.ops.div(vNW);

      if ~isempty(obj.wellModel)
        [qw,qnw] = obj.addWellFlux(mobW,mobNW,rW,rNW);
        wet  = wet  - qw;
        nwet = nwet - qnw;
      end
    end


  end


  methods (Access = protected)



  end


  methods (Access=private)


    function [wetFlux,nwetFlux] = addWellFlux(obj,mobW,mobNW,rhoW,rhoNW)

      % Signed mass sources: + injection, - production.
      % BHP is interpreted at completion elevation (no well gravity term).
      % Pure-phase injectors are supported; mixed surface compositions need PVT.
      p = obj.ops.p;
      wetFlux = 0*p;
      nwetFlux = 0*p;
      totMob = mobW + mobNW;

      for i = 1:numel(obj.wellModel.wells)
        W = obj.wellModel.wells(i);
        pbh = getValue(obj.wellModel,i,obj.getState.time);
        c = obj.domain.dofm.getLocalDoF( ...
          obj.domain.dofm.getVariableId("pressure"),W.cellList(:));
        c = c(:);
        assert(all(c > 0),"Well completions must be active pressure cells.");
        WI = W.WI(:);
        assert(isscalar(WI) || numel(WI)==numel(c), ...
          "WI must be scalar or have one entry per completion.");
        assert(all(isfinite(WI) & WI > 0),"Well indices must be positive and finite.");
        assert(isscalar(pbh) && isfinite(pbh),"Expected a finite scalar BHP.");
        phi = WI .* (pbh-p(c));

        if W.sign == 1
          comp = W.comp(:);
          assert(numel(comp)==2 && all(comp>=0) && ...
            abs(sum(comp)-1)<1e-12 && nnz(comp)==1, ...
            "Use a pure-phase injection composition, e.g. [1 0].");
          fw = comp(1).*totMob(c).*rhoW(c).*phi;
          fnw = comp(2).*totMob(c).*rhoNW(c).*phi;
        elseif W.sign == -1
          % phi is negative during production: no additional minus here.
          fw = mobW(c).*rhoW(c).*phi;
          fnw = mobNW(c).*rhoNW(c).*phi;
        else
          error("Well sign must be +1 or -1.");
        end

        % Local indexed accumulation; sum repeated completions individually.
        if numel(unique(c)) == numel(c)
          wetFlux(c) = wetFlux(c) + fw;
          nwetFlux(c) = nwetFlux(c) + fnw;
        else
          for j = 1:numel(c)
            wetFlux(c(j)) = wetFlux(c(j)) + fw(j);
            nwetFlux(c(j)) = nwetFlux(c(j)) + fnw(j);
          end
        end
      end

    end

  end





  methods (Static)
    function out = getField()
      out = ["pressure","saturation"];
      % saturation of the wetting phase
    end
  end


end

