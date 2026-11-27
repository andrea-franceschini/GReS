classdef (Abstract) Constraint < handle
  % Base class for the constraints registered on an Interface.

  properties
    name
    coupledFields
  end

  properties (GetAccess=public, SetAccess=protected)
    interface
    grid
  end

  methods (Abstract)
    registerConstraint(obj,varargin)
    assembleConstraint(obj,varargin)
    updateState(obj,solution)
    [cellData,pointData] = writeVTK(obj,interpolationFactor,t)
    writeSolution(obj,interpolationFactor,tID)
    coupledFld = getCoupledField(obj)
    fld = getField(obj)
  end

  methods
    function obj = Constraint(interface)
      if nargin > 0
        obj.interface = interface;
        obj.grid = interface.grid;
      end
    end

    function state = getState(obj,varargin)
      state = obj.interface.getState(varargin{:});
    end

    function state = getStateOld(obj,varargin)
      state = obj.interface.getStateOld(varargin{:});
    end

    function state = getStateInit(obj,varargin)
      state = obj.interface.getStateInit(varargin{:});
    end

    function setState(obj,varargin)
      obj.interface.setState(varargin{:});
    end

    function setStateOld(obj,varargin)
      obj.interface.setStateOld(varargin{:});
    end

    function setStateInit(obj,varargin)
      obj.interface.setStateInit(varargin{:});
    end

    function advanceState(obj)
      obj.setStateOld(obj.getState());
    end

    function goBackState(obj)
      obj.setState(obj.getStateOld());
    end

    function changed = updateConfiguration(obj)
      changed = false;
    end

    function resetConfiguration(obj)
    end

    function initialize(obj)
    end

    function timeStepSetup(obj)
    end

    function finalizeOutput(obj)
    end

    function out = isLinear(obj)
      out = false;
    end

    function J = getJacobian(obj,varargin)
      if isempty(varargin)
        varargin = {obj.getField()};
      end
      J = obj.interface.model.getJacobian(obj.interface,varargin{:});
    end

    function res = getResidual(obj)
      res = obj.interface.model.getResidual(obj.interface,obj.getField());
    end

    function dofm = getDoFManager(obj,varargin)
      dofm = obj.interface.getDoFManager(varargin{:});
    end

    function nDofs = getNumbDoF(obj)
      dofm = obj.getDoFManager();
      nDofs = 0;
      for field = reshape(string(obj.getField()),1,[])
        nDofs = nDofs + dofm.getNumbDoF(field);
      end
    end

    function dofs = getMultiplierDoF(obj,el)
      dofm = obj.getDoFManager();
      id = dofm.getVariableId(obj.getField());
      if nargin == 1
        dofs = dofm.getLocalDoF(id);
      else
        location = dofm.getFieldLocation(id);
        incidence = getIncidenceID(location,obj.grid,entityField.cell,el);
        dofs = dofm.getLocalDoF(id,incidence.getData());
      end
    end

    function R = getRotationMatrix(obj,side,el)
      trace = obj.grid.traceGrids(side);
      location = obj.getDoFManager().getFieldLocation(obj.getField());
      if location == entityField.cell
        R = reshape(trace.cells.rotationMatrices(el,:),3,3);
      elseif location == entityField.node
        nodes = trace.getCellNodes(el);
        normals = trace.cells.avgNodNormal(nodes,:)';
        R = zeros(numel(normals));
        for node = 1:numel(nodes)
          rows = 3*(node-1)+(1:3);
          R(rows,rows) = mxComputeRotationMat(normals(:,node));
        end
      end
    end

    function val = getVariableField(obj,side,varargin)
      field = obj.getDomainField(varargin{:});
      val = obj.interface.getVariableField(side,field);
    end

  end

  methods (Access=private)
    function field = getDomainField(obj,field)
      if nargin == 1
        field = obj.getCoupledField();
        assert(isscalar(string(field)),"Specify the coupled domain field");
      end
    end
  end
end
