classdef TwoPhaseFlow < SinglePhaseFlowFVTPFA
  % Immiscible Two-phase flow solver using MRST ADI framework
  % we assume that phases are immiscible and capillary effects are
  % negligible

  % Uses TPFA operators
  
  properties

    ops                       % struct collecting AD operators
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
      N     = G.faces.neighbors(obj.isIntFaces,:);
      cellIdx = obj.domain.dofm.getFieldCells("pressure");
      nc = numel(cellIdx);
      nf = size(N,1);

      D     = sparse([(1:nf)'; (1:nf)'], N, ones(nf,1)*[-1 1], nf, nc);
      obj.ops.grad  = @(x)  D*x;
      obj.ops.div   = @(x) -D'*x;
      obj.ops.avg   = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));
      obj.ops.upw   = @(flag, x) flag.*x(N(:, 1)) + ~flag.*x(N(:, 2));
      obj.ops.gradz = obj.ops.grad(G.cells.center(cellIdx,3));

     
    end


    function setFluidProperties(obj,input)


      defaultWet = struct('wetViscosity',[],...
        'nwetViscosity',[],...
        'refPressure',[],...
        'refWetDensity',[],...
        'refNonWetDensity',[],...
        'wetDensity',[],...
        'nwetDensity',[],...
        'wetRelPermability',[],...
        'nwetRelPermeability',[],...
        'wetCompressibility',[],...,
        'nwetCompressibility',[]);

      input = readInput(default,input);

      pR = refPressure;

      % water phase
      muW = input.wetViscosity;
      cw = input.wetCompressibility;
      rhoWR = refWetDensity;





      cw = 2.9008e-10;
      rhoWR = 1014;
      pR   = 200*barsa;
      rhoW = @(p) rhoWR .* exp( cw * (p - pR) );
      krW    = @(S) S.^2;
      obj.props.wet.density = rhoW;
      obj.props.wet.compressibility = cw;
      obj.props.wet.viscosity = muW;
      obj.props.wet.relPerm = krW;


      % oil phase
      muO = 5e-3;
      co = 1.4504e-9;
      rhoOR = 850;
      rhoO   = @(p) rhoOR .* exp( co * (p - pR) );
      krO    = @(S) S.^3;
      obj.props.nwet.density = rhoO;
      obj.props.nwet.compressibility = co;
      obj.props.nwet.viscosity = muO;
      obj.props.nwet.relPerm = krO;

      % rock properties
      dofm = obj.domain.dofm;
      regions = getTargetRegions(dofm,"pressure");
      nc = dofm.getNumbDoF("pressure");
      pv0 = zeros(nc,1);
      cr = zeros(nc,1);
      g = obj.domain.grid;
      mat = obj.domain.materials;
      pref = 2e7;


      for i = regions'
        rock = mat.getPorousRock(i);
        cIdx = g.cells.tag == i;
        poro0 = rock.getPorosity;
        pv0(cIdx) = g.cells.volume(cIdx).*poro0;
        cr(cIdx) = rock.getCompressibility;
      end

      % pore volume adi function
      obj.props.rock.poreVolume   = @(p) pv0 .* exp(cr .* (p - pref) );

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


    end



    function assembleSystem(obj,dt)

      % compute material props
      p0 = getStateOld(obj).pressure;
      sW0 = getStateOld(obj).saturation;
      p = obj.ops.p;
      sW = obj.ops.sw;

      T = obj.trans(obj.isIntFaces);


      % fluid props
      w = obj.props.wet;
      nw = obj.props.nwet;
      rW = w.density(p);
      rW0 = w.density(p0);
      rNW = nw.density(p);
      rNW0 = nw.density(p0);
      mobW = w.relPerm(sW)./w.viscosity;
      mobNW = nw.relPerm(1-sW)./nw.viscosity;

      % rock props
      pv = obj.props.rock.poreVolume;
      vol = pv(p);
      vol0 = pv(p0);

      g     = norm(gravity);
      dp  = obj.ops.grad(p);
      dpW = dp+g*obj.ops.avg(rW).*obj.ops.gradz;
      dpO = dp+g*obj.ops.avg(rNW).*obj.ops.gradz;

      % phase fluxes
      vW  = -obj.ops.upw(value(dpW) <= 0, rW.*mobW).*T.*dpW;
      vNW  = -obj.ops.upw(value(dpO) <= 0, rNW.*mobNW).*T.*dpO;
  
      % conservation of phases
      wet = (1/dt).*(vol.*rW.*sW - vol0.*rW0.*sW0) + obj.ops.div(vW);
      nwet   = (1/dt).*(vol.*rNW.*(1-sW) - vol0.*rNW0.*(1-sW0)) + obj.ops.div(vNW);


      % bcs imposition as in MRST
      % Injector: volumetric source term multiplied by surface density
      % src = 1.2860;
      % wet(1) = wet(1) - src;
      %
      % % Producer: replace equations by new ones specifying fixed pressure
      % % and zero water saturation
      % outPres = 10000000;
      % wet(end) = p(end) - outPres;
      % nwet(end)   = sW(end);

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
      p.val = p.val + dp;
      sw.val = sw.val + ds;
      % regularize saturation
      sw.val = max( min(sw.val, 1), 0);
      [obj.ops.p,obj.ops.sw] = deal(p,sw);

      setState(obj,value(p),"pressure");
      setState(obj,value(sw),"saturation");






    end



  end





  methods (Static)
    function out = getField()
      out = ["pressure","saturation"];
      % saturation of the wetting phase
    end
  end


end

