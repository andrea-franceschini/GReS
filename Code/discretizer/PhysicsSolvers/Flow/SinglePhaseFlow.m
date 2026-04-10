classdef (Abstract) SinglePhaseFlow < PhysicsSolver
  %SINGLEPHASEFLOW

  properties
    H
    P
    rhsGrav         % gravity contribution to rhs
    steadyState     % flag to force steady state problem
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

      nTags = obj.mesh.nCellTag;

      default = struct('targetRegions',1:nTags,...
        'steadyState',0);


      params = readInput(default,varargin{:});

      obj.steadyState = logical(params.steadyState);

      dofm = obj.domain.dofm;

      dofm.registerVariable(obj.getField(),fldLocation,1,params.targetRegions);
      n = getNumberOfEntities(fldLocation,obj.mesh);
      obj.fieldId = dofm.getVariableId(obj.getField());

      % initialize the state object with a pressure field
      obj.getState().data.(obj.getField()) = zeros(n,1);

    end
    

    function assembleSystem(obj,dt)

      % compute stiffness matrix obj.H and capacity matrix obj.P
      computeMat(obj,dt);

      p = obj.getState(obj.getField());
      ents = obj.domain.dofm.getActiveEntities(obj.fieldId);

      if obj.steadyState
        obj.domain.J{obj.fieldId,obj.fieldId} = obj.H;
        rhs = obj.H*p(ents);
      else
        obj.domain.J{obj.fieldId,obj.fieldId} = obj.H + obj.P/dt;
        pOld = obj.getStateOld(obj.getField());
        rhs = obj.H*p(ents) + (obj.P/dt)*(p(ents) - pOld(ents));
      end

      gamma = obj.domain.materials.getFluid().getSpecificWeight();
      if gamma > 0
        % add rhs gravity contribution
        obj.domain.rhs{obj.fieldId} = rhs + getRhsGravity(obj);
      else
        obj.domain.rhs{obj.fieldId} = rhs;  
      end
      
    end

    function updateState(obj,dSol)
      dofm = obj.domain.dofm;
      if nargin > 1
        ents = dofm.getActiveEntities(obj.fieldId);
        state = getState(obj);
        state.data.pressure(ents) = state.data.pressure(ents) + dSol(dofm.getDoF(obj.fieldId));
      end
    end

    % function advanceState(obj)
    %   % does nothing for now, but needed to override the abstract
    %   % physicsSolver method
    % end

    function [cellData,pointData] = writeVTK(obj,fac,t)
      % append state variable to output structure
      sOld = getStateOld(obj);
      sNew = getState(obj);

      p = sNew.data.pressure*fac+sOld.data.pressure*(1-fac);

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

    function alpha = getRockCompressibility(obj,el)
      mat = obj.domain.materials;
      targetRegions = getTargetRegions(obj.domain.dofm,["pressure","displacements"]);
      if ismember(obj.mesh.cellTag(el),targetRegions)
        alpha = 0; %this term is not needed in a coupled formulation
      else
        if isfield(mat.getMaterial(obj.mesh.cellTag(el)),"ConstLaw")
          %solid skeleton contribution to storage term as oedometric compressibility .
          alpha = mat.getMaterial(obj.mesh.cellTag(el)).ConstLaw.getRockCompressibility();
        else
          alpha = 0;
        end
      end
    end

    function perm = printPermeab(obj)
      % printPermeab - print the permeability for the cell or element.
      perm = zeros(obj.mesh.nCells,6);
      for el=1:obj.mesh.nCells
        ktmp = obj.domain.materials.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
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

      ssPressure = obj.domain.state.data.pressure;
      
    end
  end

  methods (Static)
    function out = getField()
      out = "pressure";
    end
  end

end
