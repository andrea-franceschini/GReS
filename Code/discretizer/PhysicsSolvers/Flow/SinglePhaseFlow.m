classdef (Abstract) SinglePhaseFlow < PhysicsSolver
  %SINGLEPHASEFLOW

  properties
    H                         % stiffness matrix
    P                         % capacity matrix
    watLev
    steadyState               % flag to force steady state problem
  end

  properties (Access = protected)
    fieldId
  end

  methods (Abstract)
    % mandatory methods that need to be implemented in any SinglePhaseFlow

    % Print the output at the physics level.
    [cellData,pointData] = buildPrintStruct(obj,outPrint);

    computeMat(obj,dt)

  end

  methods (Access = public)
    function obj = SinglePhaseFlow(domain)
      obj@PhysicsSolver(domain);
    end


    function registerSolver(obj,fldLocation,varargin)

      cells = obj.grid.cells;
      nTags = cells.nTag;

      default = struct('targetRegions',1:nTags,...
                       'steadyState',0,...
                       'waterLevel',0);


      params = readInput(default,varargin{:});

      obj.steadyState = logical(params.steadyState);
      obj.watLev = params.waterLevel;

      dofm = obj.domain.dofm;

      dofm.registerVariable(obj.getField(),fldLocation,1,params.targetRegions);
      n = getNumberOfEntities(fldLocation,obj.grid);
      obj.fieldId = dofm.getVariableId(obj.getField());

      % initialize the state object with a pressure field
      state = getState(obj);
      state.pressure = zeros(n,1);
      setState(obj,state)

    end

    function initialize(obj)
      % nothing to initialize here
    end


    function assembleSystem(obj,dt)

      % compute stiffness matrix obj.H and capacity matrix obj.P
      computeMat(obj,dt);

      p = obj.getState(obj.getField());
      ents = obj.domain.dofm.getActiveEntities(obj.fieldId);

      entType = obj.domain.dofm.getFieldLocation(obj.fieldId);
      z = getLocation(entType,obj.grid,ents);
      z = z(:,3);
      gamma = obj.domain.materials.getFluid.getSpecificWeight;

      if obj.steadyState
        obj.domain.J{obj.fieldId,obj.fieldId} = obj.H;
        rhs = obj.H*(p(ents) + gamma*z);
      else
        obj.domain.J{obj.fieldId,obj.fieldId} = obj.H + obj.P/dt;
        pOld = obj.getStateOld(obj.getField());
        rhsH = obj.H*(p(ents) + gamma*z);
        rhsP = (obj.P/dt)*(p(ents) - pOld(ents));
        rhs = rhsH + rhsP;
      end

      obj.domain.rhs{obj.fieldId} = rhs;
      
    end

    function updateState(obj,dSol)
      dofm = obj.domain.dofm;
      p = getState(obj,"pressure");
      if nargin > 1
        ents = dofm.getActiveEntities(obj.fieldId);
        p(ents) = p(ents) + dSol(dofm.getDoF(obj.fieldId));
      end
      setState(obj,p,"pressure");
    end


    function states = finalizeState(obj,p,t)
      
      % Compute the posprocessing variables for the module.
      fluid = obj.domain.materials.getFluid();
      gamma = fluid.getSpecificWeight();
      entType = obj.domain.dofm.getFieldLocation(obj.fieldId);
      z = getLocation(entType,obj.grid);
      z = z(:,3);
      mob = computeMobility(obj,p);

      if gamma>0
        states.potential = p + gamma*z;
        states.head = p/gamma + z;
      end

      states.flux = computeFlux(obj,p,mob,t);
      states.perm = printPermeab(obj);
      states.pressure = p;
    end

    % function advanceState(obj)
    %   % does nothing for now, but needed to override the abstract
    %   % physicsSolver method
    % end

    function [cellData,pointData] = writeVTK(obj,fac,t)

      p = obj.domain.state.interpolate(fac,"pressure");

      outPrint = finalizeState(obj,p,t);
      [cellData,pointData] = buildPrintStruct(obj,outPrint);
    end

    function writeSolution(obj,fac,tID)
      pOld = getStateOld(obj,obj.getField());
      pCurr = getState(obj,obj.getField());
      obj.domain.outstate.results(tID).pressure = pCurr*fac+pOld*(1-fac);
    end

    function out = isLinear(obj)
      out = true;
    end
    
    function out = isSymmetric(obj)
      out = isLinear(obj);
    end


    function alpha = getRockCompressibility(obj,cellTag,coupledRegions)

      % regions: cell tags where pressure is coupled with displacements

      mat = obj.domain.materials.getMaterial(cellTag);
      rock = mat.PorousRock;
      alpha = 0.0;

      if nargin < 3
        coupledRegions = [];
      end

      % enter only if the cell is not already coupled with displacements
      if ~ismember(cellTag,coupledRegions) 
        % get alpha from mechanical constitutive law if provided
        if isfield(mat,"ConstLaw")
          alpha = mat.ConstLaw.getRockCompressibility();
        else
          % otherwise get the compressibility from the porous rock (default
          % is 0)
          alpha = rock.getCompressibility;
        end
      end
    end

    function perm = printPermeab(obj)

      % printPermeab - print the permeability for the cell or element.
      cells = obj.grid.cells;
      perm = zeros(cells.num,6);
      for el=1:cells.num
        ktmp = obj.domain.materials.getMaterial(cells.tag(el)).PorousRock.getPermMatrix();
        perm(el,1)=ktmp(1,1);
        perm(el,2)=ktmp(2,2);
        perm(el,3)=ktmp(3,3);
        perm(el,4)=ktmp(1,2);
        perm(el,5)=ktmp(2,3);
        perm(el,6)=ktmp(1,3);
      end
    end



    function ssPressure = getSteadyStatePressure(obj,varargin)
      %
      % getSeadyStatePressure  Compute initial stress via a one-step GReS simulation.
      % DEFAULT SETUP
      % ssPressure = getSeadyStatePressure(obj)   % default setup
      % SETUP WITH USER DEFINED PARAMETERS
      % initialStress = getSeadyStatePressure(obj, 'simulationparameters', sp, 'output', out)
      %
      % Accept only a single domain with no interfaces.
      %
      % Output:
      %   initialStress - [nGP x 6] matrix with stress tensor components

      % default parameters for the simulation
      if ~obj.steadyState
        error("getSeadyStatePressure() is available only if flow solver '%s' has property 'steadyState' set to true.",class(obj))
      end

      parm = struct('Start',0.0,'End',1.0,'DtInit',1.0,'DtMax',1.0,'DtMin',0.01);
      default = struct('simulationparameters',SimulationParameters(parm),...
        'domains',obj.domain);

      params = readInput(default,varargin{:});

      solv = NonLinearImplicit(params);

      assert(solv.nDom == 1 && solv.nInterf == 0, ...
        "SteadyState pressure computation can be used only as a stand-alone PhysicsSolver within a single domain");

      solv.simulationLoop();

      ssPressure = getState(obj,obj.getField());
      
    end
  end


  methods (Access=protected)

    function [lw,dlw] = computeMobility(obj,varargin)
      mu = obj.domain.materials.getFluid().getDynViscosity();
      lw = 1/mu;
      dlw = 0.0;
    end

  end

  methods (Static)
    function out = getField()
      out = "pressure";
    end
  end

end
