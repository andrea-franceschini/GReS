classdef (Abstract) SinglePhaseFlow < PhysicsSolver
  %SINGLEPHASEFLOW

  properties
    H
    P
    rhsGrav         % gravity contribution to rhs
  end

  properties (Access = protected)
    fieldId
  end

  methods (Abstract)
    % mandatory methods that need to be implemented in any SinglePhaseFlow

    % compute the jacobian matrix
    J = computeMat(obj,dt);

    % compute the rhs
    rhs = computeRhs(obj,dt);

    % Print the output at the physics level.
    [cellData,pointData] = buildPrintStruct(obj,outPrint);
  end

  methods (Access = public)
    function obj = SinglePhaseFlow(domain)
      obj@PhysicsSolver(domain);
    end

    function assembleSystem(obj,dt)
      obj.domain.J{obj.fieldId,obj.fieldId} = computeMat(obj,dt);
      obj.domain.rhs{obj.fieldId} = computeRhs(obj,dt);
    end

    function updateState(obj,dSol)
      if nargin > 1
        ents = obj.dofm.getActiveEntities(obj.fieldId);
        state = getState(obj);
        state.data.pressure(ents) = state.data.pressure(ents) + dSol(obj.dofm.getDoF(obj.fieldId));
      end
    end

    function advanceState(obj)
      % does nothing for now, but needed to override the abstract
      % physicsSolver method
    end

    function [cellData,pointData] = writeVTK(obj,t)
      % append state variable to output structure
      sOld = getStateOld(obj);
      sNew = getState(obj);

      if isempty(sOld)
        p = sNew.data.pressure;
      else
        fac = (t - sOld.t)/(sNew.t - sOld.t);
        p = sNew.data.pressure*fac+sOld.data.pressure*(1-fac);
      end

      outPrint = finalizeState(obj,p,t);
      [cellData,pointData] = buildPrintStruct(obj,outPrint);
    end

    function writeMatFile(obj,t,tID)
      pOld = getStateOld(obj,obj.getField());
      pCurr = getState(obj,obj.getField());
      fac = (t - getStateOld(obj).t)/(getState(obj).t - getStateOld(obj).t);
      obj.domain.outstate.results(tID).pressure = pCurr*fac+pOld*(1-fac);
    end

    function out = isLinear(obj)
      out = true;
    end

    function alpha = getRockCompressibility(obj,el)
      if ismember(obj.mesh.cellTag(el),getTargetRegions(obj.dofm,["pressure","displacements"]))
        alpha = 0; %this term is not needed in a coupled formulation
      else
        if isfield(obj.materials.getMaterial(obj.mesh.cellTag(el)),"ConstLaw")
          %solid skeleton contribution to storage term as oedometric compressibility .
          alpha = obj.materials.getMaterial(obj.mesh.cellTag(el)).ConstLaw.getRockCompressibility();
        else
          alpha = 0;
        end
      end
    end

    function perm = printPermeab(obj)
      % printPermeab - print the permeability for the cell or element.
      perm = zeros(obj.mesh.nCells,6);
      for el=1:obj.mesh.nCells
        ktmp = obj.materials.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
        perm(el,1)=ktmp(1,1);
        perm(el,2)=ktmp(2,2);
        perm(el,3)=ktmp(3,3);
        perm(el,4)=ktmp(1,2);
        perm(el,5)=ktmp(2,3);
        perm(el,6)=ktmp(1,3);
      end
    end
  end

  methods (Static)
    function out = getField()
      out = "pressure";
    end
  end

end
