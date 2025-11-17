classdef (Abstract) physicsSolver < handle
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
  %   J           - System Jacobian (abstract, must be implemented by
  %   subclass). rhs         - Right-hand side vector (abstract, must be
  %   implemented by subclass).
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


  properties (Abstract, GetAccess=public, SetAccess=private)
    J
    rhs
  end

  properties (GetAccess=public, SetAccess=protected)
    fields
  end


  properties (GetAccess=public, SetAccess=protected)
    dofm                % handle to dofManager
    model               % handle to entire model
    bcs                 % boundary conditions
    outstate            % printing utilities
    materials           % materials
    mesh                % topology of the grid
    elements            % db for finite elements
    faces               % db for finite volumes
    variables           % list of fields with the same order of the global system
  end

  properties
    state               % the state object holding all variables in the model
    stateOld            % the state object after the last converged time step
  end

  methods
    function obj = physicsSolver(model,grid,mat,bcs,inputStruct)

      % inputStruct: struct with additional solver-specific parameters

      obj.model = model;
      obj.mesh = grid.topology;
      obj.elements = grid.elements;
      obj.faces = grid.faces;
      obj.materials = mat;
      obj.bcs = bcs;

      % create the DoFManager
      obj.dofm = DoFManagerNew(grid);

      % create the State object
      obj.state = State();

      % read xml fields that must have the same name of the class itself
      % this is where the reading solver specific input from file
      obj.registerSolver(solverInput);

    end
  end

  methods (Abstract)

    % mandatory methods that need to be implemented in the physicsSolver

    % read the input data of the solver and assign variables to cell tags
    registerSolver(obj,input);

    % compute the jacobian and the rhs
    assembleSystem(obj);

    % get the boundary condition dofs and vals (solver specific)
    [bcDofs,bcVals] = getBC(obj,bcId);

    % update the state variables after solving the problem
    updateState(obj,solution);

    % advance the state variables after achieving solver convergence
    advanceState(obj);

    % update the output structres for printing purposes
    [cellData,pointData] = printState(obj,t);
  end


  methods

    function J = getJacobian(obj,varargin)
      % GETJACOBIAN Return the system Jacobian matrix
      %
      % Usage:
      %   J = getJacobian(obj)
      %       Returns the full Jacobian matrix.
      %
      %   J = getJacobian(obj, fieldList)
      %       Returns only the Jacobian blocks corresponding to the
      %       specified fields. Assumes the same fields for both rows and columns.
      %
      %   J = getJacobian(obj, rowFields, colFields)
      %       Returns the Jacobian blocks corresponding to the specified
      %       fields for rows and columns separately.
      %
      % Inputs:
      %   fieldList  - string or cell array of field names (for both rows and columns)
      %   rowFields  - string or cell array of field names for rows
      %   colFields  - string or cell array of field names for columns
      %
      % Output:
      %   J          - the assembled Jacobian matrix
      %
      % Notes:
      %   - If no input fields are provided, the full Jacobian is returned.
      %   - Row and column fields must correspond to existing variables in the system.

      if nargin == 1
        J = obj.J;
      elseif nargin == 2
        id = obj.dofm.getVariableId(varargin{1});
        dofList = getDoF(obj.dofm,id);
        J = obj.J(dofList,dofList);
      elseif nargin == 3
        idRow = obj.dofm.getVariableId(varargin{1});
        idCol = obj.dofm.getVariableId(varargin{2});
        dofsRow = getDoF(obj.dofm,idRow);
        dofsCol = getDoF(obj.dofm,idCol);
        obj.dofm.getVariableId(varargin{2});
        J = obj.J(dofsRow,dofsCol);
      else
        error("Too many input arguments")
      end

    end

    function rhs = getRhs(obj,varargin)
      % GETRHS Return the right-hand side vector
      %
      % Usage:
      %   rhs = getRhs(obj)               - returns full RHS
      %   rhs = getRhs(obj, fieldList)    - returns only specified fields
      %
      % Inputs:
      %   fieldList - string or cell array of field names
      %
      % Output:
      %   rhs       - RHS vector (subset or full)
      %
      % Notes:
      %   Only one field list is allowed; multiple fields will be concatenated

      if nargin == 1
        rhs = obj.rhs;
      elseif nargin == 2
        id = obj.dofm.getVariableId(varargin{1});
        dofList = getDoF(obj.dofm,id);
        rhs = obj.rhs(dofList);
      else
        error("Too many input arguments")
      end
    end

 
    function applyDirVal(obj,bcDofs,bcVals,bcVariableName)

      % apply Dirichlet BC values to state variables
      obj.state.data.(bcVariableName)(bcDofs) = bcVals;
    end

    function applyBCs(obj,t)
      bcList = keys(obj.bcs.keys);

      for bcId = string(bcList)
        [bcDofs,bcVals] = getBC(obj,bcId,t); 
        % TO DO: use model property to remove slave dofs in mortar
        bcVar = obj.bcs.getVariable(bcId);
        obj.applyBC(bcId,bcDofs,bcVals,bcVar)
      end
    end

    function applyBC(obj,bcId,bcDofs,bcVals,bcVar)

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

      obj.rhs{bcId}(bcDofs) = obj.rhs{bcId}(bcDofs) - bcVals;

    end

    function applyDirBC(obj,bcDofs,bcVariableName)

      % Standard application of Dirichlet boundary condition to the jacobian.
      % This method works with incremental linear system du = J\(-rhs)

      % sort bcDofs to improve sparse access performance
      bcDofs = sort(bcDofs);

      bcId = obj.dofm.getVariableId(bcVariableName);

      % zero out columns
      for i = 1:numel(fields)
        obj.J{i,bcId}(:,bcDofs) = 0;
      end

      % zero out rows (use transpose trick)
      for j = 1:numel(fields)
        obj.J{bcId,j} = obj.J{bcId,j}';
        obj.J{bcId,j}(:,bcDofs) = 0;
        obj.J{bcId,j} = obj.J{bcId,j}';
      end

      % add 1 to diagonal entry of diagonal block
      obj.J{bcId,bcId}(sub2ind(size(obj.J{{bcId,bcId}}), bcDofs, bcDofs)) = 1;

    end
  end
end
