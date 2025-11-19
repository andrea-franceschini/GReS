classdef (Abstract) PhysicsSolver < handle
  % PHYSICSSOLVER Abstract interface for creating physics solvers in GReS
  %
  % This abstract class defines the interface for any physics solver in
  % GReS. Any subclass must implement all functionality required to
  % assemble the full Jacobian and the right-hand side (RHS) vector of the
  % simulation.
  %
  % Properties:
  %   fields      - String array of variable fields that are going to be used by
  %   the solver
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
    dofm
    mesh
    elements
    faces
    bcs
    materials
    simparams
  end

  methods
    function obj = PhysicsSolver(domain)

      % inputStruct: struct with additional solver-specific parameters
      obj.domain = domain;
      obj.dofm = domain.dofm;
      obj.mesh = domain.grid.topology;
      obj.elements = domain.grid.cells;
      obj.faces = domain.grid.faces;
      obj.materials = domain.materials;
      obj.bcs = domain.bcs;
      obj.simparams = domain.simparams;

    end
  end

  methods (Abstract)

    % mandatory methods that need to be implemented in the physicsSolver

    % read the input data of the solver and assign variables to cell tags
    registerSolver(obj,input);

    % compute the jacobian and the rhs
    assembleSystem(obj,varargin);

    % get the boundary condition dofs and vals (solver specific)
    %[bcDofs,bcVals] = getBC(obj,bcId,t);

    % update the state variables after solving the linear system
    updateState(obj,solution);

    % update the state variables after achieving solver convergence
    advanceState(obj);

    % update the output structres for printing purposes
    [cellData,pointData] = printState(obj,t);
  end


  methods

    % interface to get and set the state object from the solver

    function stat = getState(obj,varName)
      % get a copy of a state variable field
      if nargin < 2
        stat = obj.domain.getState();
      else
        stat = obj.domain.getState().data.(varName);
      end
    end

    function stat = getStateOld(obj,varName)
      % get a copy of a state variable field
      if isempty(varName)
        stat = obj.domain.getStateOld();
      else
        stat = obj.domain.getStateOld().data.(varName);
      end
    end

 
    function applyDirVal(obj,bcDofs,bcVals,bcVariableName)

      % apply Dirichlet BC values to state variables
      obj.setState(bcVariableName,bcDofs,bcVals);

    end

    function applyBCs(obj,t)
      bcList = keys(obj.bcs.keys);

      for bcId = string(bcList)
        bcVar = obj.bcs.getVariable(bcId);

        if any(strcmpi(obj.fields,bcVar))

          % TO DO: use model property to remove slave dofs in mortar
          obj.applyBC(t,bcId,bcVar)
        end
      end
    end

    % implement here a base getBC method
    % function [bcDofs,bcVals] = getBC(obj,bcId,t)
    % 
    %   error("Error in %s: Boundary condition type '%s' is not " + ...
    %     "available in the base version of applyBC(). Consider " + ...
    %     "overriding this method for a specific implementation", ...
    %     class(obj),bcType);
    % 
    % end

    function applyBC(obj,t,bcId,bcVar)

      % get bcDofs and bcVals (should be implemented by the physicsSolver)
      [bcDofs,bcVals] = getBC(obj,bcId,t);

      % Base application of a Boundary condition
      bcType = obj.bcs.getType(bcId);

      switch bcType
        case 'Dirichlet'
          applyDirBC(obj,bcDofs,bcVals,bcVar);
        case 'Neumann'
          applyNeuBC(obj,bcDofs,bcVar);
        otherwise
          error("Error in %s: Boundary condition type '%s' is not " + ...
            "available in the base version of applyBC(). Consider " + ...
            "overriding this method for a specific implementation", ...
            class(obj),bcType);
      end
    end

    function applyNeuBC(obj,bcDofs,bcVals,bcVariableName)

      % Base application of Neumann boundary condition to the rhs.
      % bc values are subtracted since we solve du = J\(-rhs)
      bcId = obj.dofm.getVariableId(bcVariableName);

      obj.domain.rhs{bcId}(bcDofs) = obj.domain.rhs{bcId}(bcDofs) - bcVals;

    end

    function applyDirBC(obj,bcDofs,bcVariableName)

      % Standard application of Dirichlet boundary condition to the jacobian.
      % This method works with incremental linear system du = J\(-rhs)

      % sort bcDofs to improve sparse access performance
      bcDofs = sort(bcDofs);

      bcId = obj.dofm.getVariableId(bcVariableName);

      % zero out columns
      nV = getNumberOfVariables(obj.dofm);

      for i = 1:nV
        obj.domain.J{i,bcId}(:,bcDofs) = 0;
      end

      % zero out rows (use transpose trick)
      for j = 1:nV
        obj.domain.J{bcId,j} = obj.J{bcId,j}';
        obj.domain.J{bcId,j}(:,bcDofs) = 0;
        obj.domain.J{bcId,j} = obj.J{bcId,j}';
      end

      % add 1 to diagonal entry of diagonal block
      obj.domain.J{bcId,bcId}...
        (sub2ind(size(obj.domain.J{{bcId,bcId}}), bcDofs, bcDofs)) = 1;

      % rhs treatment

    end
  end
end
