classdef (Abstract) PhysicsSolver < handle
  % PHYSICSSOLVER Abstract interface for creating physics solvers in GReS
  %
  % This abstract class defines the interface for any physics solver in
  % GReS. Any subclass must implement all functionality required to
  % assemble the full Jacobian and the right-hand side (RHS) vector of the
  % simulation.
  %
  % Properties:
  %   fields      - String array of variable fields that are going to be
  %   used by the solver
  %
  %
  % Coupled Solver Notes:
  %   - If the solver is coupled, 'fields' lists all the single physics
  %   fields.
  % - 'J' and 'rhs' are cell arrays, each cell corresponding to a
  %   single field according to the DoFManager order.
  %
  %
  % Constructor:
  %   obj = physicsSolver(inputStruct)
  %       Reads input configuration from XML. The input must have a field
  %       with the same name as the class itself.


  properties (Abstract)
    % the fields modified by the solver
    %fields
  end


  properties (GetAccess=public, SetAccess=protected)
    % handle to domain properties
    domain
    simparams
    grid
  end

  methods
    function obj = PhysicsSolver(domain)

      % inputStruct: struct with additional solver-specific parameters
      obj.domain = domain;
      obj.simparams = domain.simparams;
      obj.grid = domain.grid;

    end
  end

  methods (Abstract)
    % mandatory methods that need to be implemented in any physicsSolver

    % read the input data of the solver and assign variables to cell tags
    registerSolver(obj,varargin);

    % compute the jacobian and the rhs
    assembleSystem(obj,varargin);

    % apply the boundary condition to the jacobian and rhs
    %applyBC(obj,bcId,t);

    % apply the dirichlet values to the state object
    %applyDirVal(obj,bcId,t);
    
    % update the state variables after solving the linear system
    updateState(obj,solution);

    % update the output structures for printing purposes
    [cellData,pointData] = writeVTK(obj,interpolationFactor,t);

    % write history to MAT-file
    writeSolution(obj,interpolationFactor,tID);
  end

  methods (Abstract, Static)
    % get the list of variable fields affected by the solver

    getField();
  end


  methods
    % interface to get and set the state object from the solver

    function state = getState(obj,varargin)
      % get a copy of a state variable field
      state = obj.domain.getState(varargin{:});
    end

    function state = getStateOld(obj,varargin)
      % get a copy of a state variable field
      state = obj.domain.getStateOld(varargin{:});
    end

    function state = getStateInit(obj,varargin)
      % get a copy of a state variable field
      state = obj.domain.getStateInit(varargin{:});
    end

    function setState(obj,varargin)
      % get a copy of a state variable field
      obj.domain.setState(varargin{:});
    end

    function setStateOld(obj,varargin)
      % get a copy of a state variable field
      obj.domain.setStateOld(varargin{:});
    end

    function setStateInit(obj,varargin)
      % get a copy of a state variable field
      obj.domain.setStateInit(varargin{:});
    end

    function advanceState(obj,varargin)
      % base method to advance the state after reaching convergence
      % copy variables from current state to old state
      setStateOld(obj,getState(obj));
    end

    function goBackState(obj,varargin)
      % base method to move back the state when convergence is not reached

      setState(obj,getStateOld(obj));
    end

    function hasConfigurationChanged = updateConfiguration(obj)
      % base physicsSolver class implements no configuration update
      hasConfigurationChanged = false;
    end

    function resetConfiguration(obj)
      % reset the physicsSolver configuraiton
    end

    function initialize(obj)
      % initialize the physics solver just before the simulation starts
    end

    function timeStepSetup(obj)
      % prepare the physics solver at the begin of each time step
    end

    function finalizeOutput(obj)
      % override this method in a PhysicsSolver to produce additional
      % output files other than the general Discretizer.outState pvd file
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
    end

    function applyBC(obj,bcId,varargin)
      % applyBC apply a boundary condition to the system
      %
      %   applyBC(obj, bcId, t) applies the boundary condition identified
      %   by bcId and retreive bc dofs and values at time t
      %
      %   applyBC(obj, bcId, dofs, vals) applies the boundary condition
      %   using explicitly provided degrees of freedom dofs and
      %   corresponding values vals.
      %
      %   - for custom boundary conditions, dofs and values must be
      %   provided explicitely

      bcs = obj.domain.bcs;

      if ~BCapplies(obj,bcId)
        return
      end

      if nargin == 3
        type = bcs.getType(bcId);
        assert(~BCtype.isCustomBC(type),"Base applyBC method for custom BC requires " + ...
          "dofs and values to be already specified")
        t = varargin{1};
        dofs = bcs.getDofs(bcId,obj.domain.dofm);
        vals = bcs.getVals(bcId,t);
      elseif nargin==4
        dofs = varargin{1};
        vals = varargin{2};
      else
        error("Wrong number of input for applyBC function")
      end

      % standard application of a boundary condition

      if isEssential(bcs,bcId)
        obj.applyDirBC(bcId,dofs,vals);
      else
        obj.applyNeuBC(bcId,dofs,vals);
      end
    end

    function applyNeuBC(obj,bcId,bcDofs,bcVals)
      if ~BCapplies(obj,bcId)
        return
      end

      domain = obj.domain;
      field = domain.bcs.getVariable(bcId);
      dofs = domain.model.getGlobalDoFs(domain,field);
      active = bcDofs ~= 0;
      domain.model.addResidualAtDoFs(dofs(bcDofs(active)),-bcVals(active));

    end


    function applyDirBC(obj,bcId,bcDofs,varargin)

      if ~BCapplies(obj,bcId)
        return
      end

      domain = obj.domain;
      field = domain.bcs.getVariable(bcId);
      domain.model.applyDirichlet(domain,field,bcDofs(bcDofs ~= 0),obj.isSymmetric());

    end


    function J = getJacobian(obj,varargin)

      if isempty(varargin)
        varargin = {obj.getField()};
      end
      J = obj.domain.model.getJacobian(obj.domain,varargin{:});

    end


    function out = BCapplies(obj,bcId)
      bcVar = obj.domain.bcs.getVariable(bcId);
      out = any(strcmp(obj.getField(),bcVar));
    end

    function out = isSymmetric(obj)
      out = false;
      % optional solver query to know if a solver is symmetric or not
    end

    function out = isLinear(obj)
      out = false;
      % optional solver query to know if a solver is symmetric or not
    end

  end


  methods (Static)

  end

end
