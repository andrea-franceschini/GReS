classdef (Abstract) physicsSolver < handle
  % PHYSICSSOLVER Abstract interface for creating physics solvers in GReS
  %
  % This abstract class defines the interface for any physics solver in
  % GReS. Any subclass must implement all functionality required to
  % assemble the full Jacobian and the right-hand side (RHS) vector of the
  % simulation.
  %
  % Properties:
  %   fields      - List of single physics modules implemented in the
  %   solver.
  %                 For coupled solvers, this is a string array containing
  %                 all individual physical fields.
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
    dofm        % handle to dofManager
    simparams   % handle to simulation parameters
    bcs         % boundary conditions
    outstate    % printing utilities
    materials   % materials
    grid        % topology, elements, faces
    variables   % list of fields with the same order of the global system
  end

  methods
    function obj = physicsSolver(inputStruct)

      % domain:  handle to the Discretizer object storing all the
      % information of the model


      % read xml field that must have the same name of the class itself
      obj.readInput(inputStruct.(class(obj)));
    end
  end

  methods (Abstract)

    readInput(obj,input);
    setup(obj); 
    % register the existing state variables and DoF fields
    assembleSystem(obj); % compute the matrix and the rhs
    updateState(obj); % update the state variables 
    finalizeState(obj); % for post a processing variable
    [cellData,pointData] = buildPrintStruct(obj); % save outuputs to structures

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

    function applyDirBC(obj,bcDofs,bcVariable,colVariable)
      % Base application of boundary condition to the jacobian. This
      % method only implements standard Dirichlet BCs on diagonal jacobian
      % blocks, but can be overridden by specific implementations. This
      % version works with incremental linear system (vals = 0). BC values
      % for Dirichlet BC are zero (since the system is solved in
      % incremental form) set Dir rows to zero

      idRow = obj.dofm.getVariableId(bcVariable);

      if isempty(colVariable)
        idCol = idRow;
      else
        idCol = obj.dofm.getVariableId(colVariable);
      end


      
      obj.J{idRow,idCol} = obj.J{idRow,idCol}';
      obj.J{idRow,idCol}(:,bcDofs) = 0; % setting columns is much faster
      obj.J{idRow,idCol} = obj.J{idRow,idCol}';
      % Update rhs with columns to be removed
      %obj.rhs = obj.rhs - obj.J(:,dofs)*vals;
      % set obj.rhs vals to dir vals
      obj.rhs(bcDofs) = 0;
      obj.J{idRow,idCol}(:,bcDofs) = 0;
      % modify diagonal entries of K (avoiding for loop and accessing)
      Jdiag = diag(obj.J{idRow,idCol});
      Jdiag(bcDofs) = 1;
      obj.J{idRow,idCol} = obj.J{idRow,idCol} - diag(diag(obj.J{idRow,idCol})) + diag(Jdiag);
    end

  end
end
