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

    % mandatory methods that need to be implemented in any physicsSolver

    % read the input data of the solver and assign variables to cell tags
    registerSolver(obj,input);

    % compute the jacobian and the rhs
    assembleSystem(obj,varargin);

    % apply the boundary condition to the jacobian and rhs
    applyBC(obj,bcId,t);

    % apply the dirichlet values to the state object
    applyDirVal(obj,bcId,t);
    
    % update the state variables after solving the linear system
    updateState(obj,solution);

    % update the output structures for printing purposes
    [cellData,pointData] = writeVTK(obj,interpolationFactor,t);

    % write history to MAT-file
    writeMatFile(obj,interpolationFactor,tID);

  end

  methods (Abstract, Static)

    % get the list of variable fields affected by the solver
    getField();

  end


  methods

    % interface to get and set the state object from the solver

    function stat = getState(obj,varName)
      % get a copy of a state variable field
      if nargin < 2
        stat = obj.domain.getState();
      else
        if ~isfield(obj.domain.getState().data,varName)
          error("Variable %s does not exist in the State object",varName)
        end
        stat = obj.domain.getState().data.(varName);
      end
    end

    function stat = getStateOld(obj,varName)
      % get a copy of a state variable field
      if nargin < 2
        stat = obj.domain.getStateOld();
      else
        if ~isfield(obj.domain.getStateOld().data,varName)
          error("Variable %s does not exist in the StateOld object",varName)
        end
        stat = obj.domain.getStateOld().data.(varName);
      end
    end

    function advanceState(obj,varargin)

      % base method to advance the state after reaching convergence
      % hard copy the new state object

      obj.domain.stateOld = copy(obj.domain.state);

    end

    function goBackState(obj,varargin)
      % base method to move back the state when convergence is not reached

      obj.domain.state = copy(obj.domain.stateOld);

    end


    function applyNeuBC(obj,bcId,bcDofs,bcVals)

      if ~BCapplies(obj,bcId)
        return
      end

      % Base application of Neumann boundary condition to the rhs.
      % bc values are subtracted since we solve du = J\(-rhs)
      bcVar = obj.bcs.getVariable(bcId);
      bcId = obj.dofm.getVariableId(bcVar);

      obj.domain.rhs{bcId}(bcDofs) = obj.domain.rhs{bcId}(bcDofs) - bcVals;

    end

    function applyDirBC(obj,bcId,bcDofs,bcVals)

      % Standard application of Dirichlet boundary condition to the jacobian.
      % This method works with incremental linear system du = J\(-rhs)

      if ~BCapplies(obj,bcId)
        return
      end

      % sort bcDofs to improve sparse access performance
      [bcDofs,sortId] = sort(bcDofs);

      bcVar = obj.bcs.getVariable(bcId);
      bcVarId = obj.dofm.getVariableId(bcVar);

      nV = getNumberOfVariables(obj.dofm);

      % zero out rows (use transpose trick)
      for j = 1:nV
        obj.domain.J{bcVarId,j} = obj.domain.J{bcVarId,j}';
        obj.domain.J{bcVarId,j}(:,bcDofs) = 0;
        obj.domain.J{bcVarId,j} = obj.domain.J{bcVarId,j}';
      end

      % apply BC to multi-domain jacobian coupling blocks
      for iI = 1:numel(obj.domain.interfaces)
        if ~isempty(obj.domain.Jum{bcVarId})
          obj.domain.Jum{iI}{bcVarId} = obj.domain.Jum{iI}{bcVarId}';
          obj.domain.Jum{iI}{bcVarId}(:,bcDofs) = 0;
          obj.domain.Jum{iI}{bcVarId} = obj.domain.Jum{iI}{bcVarId}';
        end
      end

      % zero out columns only if the solver is symmetric (preserves
      % symmetry)
      % if isSymmetric(obj)
      %   for i = 1:nV
      %   obj.domain.J{i,bcVarId}(:,bcDofs) = 0;
      % end
      % 
      % for iI = 1:numel(obj.domain.interfaces)
      %   if ~isempty(obj.domain.Jum{bcVarId})
      %     obj.domain.Jmu{iI}{bcVarId}(:,bcDofs) = 0;
      %   end
      % end

      % add 1 to diagonal entry of diagonal block
      %J(bcDofs + (bcDofs-1)*size(J,1)) = 1;   extremely slow

      J = obj.domain.J{bcVarId,bcVarId};
      J = J + sparse(bcDofs, bcDofs, ones(size(bcDofs)), size(J,1), size(J,2));
      obj.domain.J{bcVarId,bcVarId} = J;

      if nargin > 3
        obj.domain.rhs{bcVarId}(bcDofs) = bcVals(sortId);
      else
        obj.domain.rhs{bcVarId}(bcDofs) = 0;
      end

    end

    function J = getJacobian(obj,varargin)

      % get the Jacobian blocks affected by the solver
      % differently from the getJacobian() in Discretizer, this method
      % returns a matrix ready to perform computations

      if nargin < 2
        J = obj.domain.getJacobian(obj.getField());
      else
        J = obj.domain.getJacobian(varargin{:});
      end

      J = cell2matrix(J);

    end

    function out = BCapplies(obj,bcId)

      bcVar = obj.bcs.getVariable(bcId);
      out = any(strcmp(obj.getField(),bcVar));

    end
    

  end


  methods (Static)

    function out = isSymmetric()

      out = false;
      % optional solver query to know if a solver is symmetric or not

    end
  end

end
