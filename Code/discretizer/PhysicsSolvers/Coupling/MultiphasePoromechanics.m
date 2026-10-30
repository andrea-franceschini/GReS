classdef MultiphasePoromechanics < BiotFullyCoupled
  % Biot model subclass
  % Coupled Poromechanics with SinglePhaseFlow
  % work only with unitary biot coefficient


  properties (Access = protected)

    % we avoid multiple inheritance and we directly create instances to the
    % single physics models that are needed


    fldP      % wetting phase equation (pressure)
    fldS      % non wetting phase equation (saturation)

  end

  methods (Access = public)
    function obj = MultiphasePoromechanics(domain)

      obj@BiotFullyCoupled(domain);

    end

    function registerSolver(obj,varargin)

      registerSolver@BiotFullyCoupled(obj,'flowSolver',"TwoPhaseFlow",varargin{:});

      % hellper
      dofm = obj.domain.dofm;
      obj.fldP =  dofm.getVariableId("pressure");
      obj.fldS =  dofm.getVariableId("saturation");


    end

    function assembleSystem(obj,dt,varargin)

      % get the mechanical Auu block
      obj.mechSolver.assembleSystem(dt);

      % biot coefficient is assumed constant = 1
      computeMat(obj);

      % assign coupling blocks to jacobian
      obj.domain.J{obj.fldMech,obj.fldP} = -obj.Q;


      dofm = obj.domain.dofm;
      flowEnts = dofm.getFieldCells("pressure");
      mechEnts = dofm.getActiveEntities("displacements",1);

      curr  = getState(obj);
      prev = getStateOld(obj);

      uOld = prev.displacements(mechEnts);
      pOld  = prev.pressure(flowEnts);
      sWOld = prev.saturation(flowEnts);

      p0 = getStateInit(obj).pressure;
     

      [u,p,sW] = initVariablesADI( ...
        curr.displacements(mechEnts), ...
        curr.pressure(flowEnts), ...
        curr.saturation(flowEnts));
      obj.flowSolver.ops.u  = u;

      rhsMech = -obj.Q * (curr.pressure(flowEnts) - p0(flowEnts));

      pv = obj.flowSolver.props.rock.poreVolume;

      [wet,nwet] = obj.flowSolver.computeFlowResidual( ...
        dt,p,sW,pOld,sWOld,pv(u,p),pv(uOld,pOld));


      % residual equations
      eqIdx = [dofm.getVariableId("displacements"); ...
               dofm.getVariableId("pressure"); ...
               dofm.getVariableId("saturation")];

      obj.domain.rhs{eqIdx(1)} = obj.domain.rhs{obj.fldMech} + rhsMech;
      obj.domain.rhs{eqIdx(2)} = wet.val;
      obj.domain.rhs{eqIdx(3)} = nwet.val;
      obj.domain.J(eqIdx(2),eqIdx) = wet.jac;
      obj.domain.J(eqIdx(3),eqIdx) = nwet.jac;

    end

    function initialize(obj)

      initialize@BiotFullyCoupled(obj);


      computeMat(obj);


      obj.flowSolver.props.rock.poreVolume = getPoreVolume(obj);

    end

    function computeMat(obj)

      % call method according to the discretization technique chosen
      if isempty(obj.Q) || ~isLinear(obj)
        computeMatBiot(obj)
      end
    end


    function out = isLinear(obj)
      out = false;
    end

    function out = getFlowScheme(obj)
      out = obj.flowSolver.typeDiscretization;
    end


    function updateState(obj,dsol)


      updateState@BiotFullyCoupled(obj,dsol);

      % update the adi variable as well as the state object
      dofm = obj.domain.dofm;
      du = dsol(dofm.getDoF(dofm.getVariableId("displacements")));
      u = obj.flowSolver.ops.u;
      u.val = u.val + du;
      obj.flowSolver.ops.u = u;
  
    end
  end


  methods (Access=protected)

    function pV = getPoreVolume(obj)

      % return function handle with pore volume expression

      % rock properties
      dofm = obj.domain.dofm;
      mechEnts = dofm.getActiveEntities("displacements",1);
      u0 = obj.getStateInit.displacements(mechEnts);
      g = obj.domain.grid;
      regions = getTargetRegions(dofm,"pressure");
      nc = dofm.getNumbDoF("pressure");
      pv0 = zeros(nc,1);
      mat = obj.domain.materials;

      for i = regions'
        rock = mat.getPorousRock(i);
        cIdx = g.cells.tag == i;
        poro0 = rock.getPorosity;
        pv0(cIdx) = g.cells.volume(cIdx).*poro0;
      end

      % 0 * p because Biot coefficient is assumed unitary
      pV = @(u,p) pv0(cIdx) + obj.Q' * (u - u0) + 0 * p .* g.cells.volume(dofm.getFieldCells("pressure"));

    end

  end

  methods (Static)

    function out = getField()
      out = [Poromechanics.getField(), TwoPhaseFlow.getField()];
    end

    function out = isSymmetric()
      out = false;
    end

  end

end


