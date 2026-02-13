classdef Discretizer < handle
  % General discretizer class

  properties (GetAccess=public, SetAccess=public)
    physicsSolvers                     % physics solvers database
    solverNames
    dofm
    bcs
    materials
    grid
  end

  properties
    simparams
    J
    rhs
    outstate
    domainId
    % Jacobian blocks for multidomain coupling
    Jum
    Jmu
    vtmBlock
  end

  properties (GetAccess=public, SetAccess=public)

    % handle to interface solver defined on the domain
    interfaces

    % id of the interfaces defined on the domain
    interfaceList

    state
    stateOld
  end


  methods (Access = public)
    function obj = Discretizer(varargin)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here

      % Initialized internal variables
      obj.dofm = DoFManager();
      obj.bcs = Boundaries();
      obj.outstate = OutState();
      obj.materials = Materials();
      obj.grid = struct('topology',Mesh(),'faces',[],'cells',[]);

      obj.physicsSolvers = containers.Map('KeyType','char','ValueType','any');
      obj.setDiscretizer(varargin{:});
    end

    function applyBC(obj,t)
      bcList = obj.bcs.getBCList();

      for bcId = bcList
        % loop over available bcs
        bcVar = obj.bcs.getVariable(bcId);

        for solv = obj.solverNames
          % loop over available solvers
          solver = obj.getPhysicsSolver(solv);

          % check if BCs applies to the solver
          if any(strcmpi(solver.getField(),bcVar))
            solver.applyBC(bcId,t)
          end

        end
      end

    end


    function applyDirVal(obj,t)
      bcList = obj.bcs.getBCList;

      for bcId = bcList
        % discard non-Dirichlet BC
        if ~strcmp(obj.bcs.getType(bcId),"Dirichlet")
          continue
        end

        bcVar = obj.bcs.getVariable(bcId);

        for solv = obj.solverNames
          % loop over available solvers
          solver = obj.getPhysicsSolver(solv);

          % check if BCs applies to the solver
          if any(strcmpi(solver.getField(),bcVar))
            solver.applyDirVal(bcId,t)
          end
        end
      end
    end

    function updateState(obj,du)

      for solv = obj.solverNames
        % loop over available solvers
        obj.getPhysicsSolver(solv).updateState(du);
      end

    end

    function isConfigChanged = updateConfiguration(obj)

      isConfigChanged = false;

      for solv = obj.solverNames
        % loop over available solvers
        isConfigChanged = any([isConfigChanged; ...
          obj.getPhysicsSolver(solv).updateConfiguration()]);
      end

    end

    function resetConfiguration(obj)


      for solv = obj.solverNames
        % loop over available solvers
        obj.getPhysicsSolver(solv).resetConfiguration();
      end

    end


    function advanceState(obj)

      for solv = obj.solverNames
        % loop over available solvers
        obj.getPhysicsSolver(solv).advanceState();
      end

    end

    function goBackState(obj)

      for solv = obj.solverNames
        % loop over available solvers
        obj.getPhysicsSolver(solv).goBackState();
      end

    end

    function finalizeOutput(obj)

      obj.outstate.finalize();

      for solv = obj.solverNames
        % loop over available solvers
        obj.getPhysicsSolver(solv).finalizeOutput();
      end

    end


    function out = getPhysicsSolver(obj,id)
      out = obj.physicsSolvers(id);
    end

    function stat = getState(obj,varName)
      % get a copy of a state variable field
      if nargin < 2
        stat = obj.state;
      else
        if ~isfield(obj.getState().data,varName)
          error("Variable %s does not exist in the State object",varName)
        end
        stat = obj.getState().data.(varName);
      end
    end

    function stat = getStateOld(obj,varName)
      % get a copy of a state variable field
      if nargin < 2
        stat = obj.stateOld;
      else
        if ~isfield(obj.getStateOld().data,varName)
          error("Variable %s does not exist in the StateOld object",varName)
        end
        stat = obj.getStateOld().data.(varName);
      end
    end


    % function addInterface(obj,interfId,interf)
    %   % add mortar interface to current domain
    %   if ~ismember(interfId,obj.interfaceList)
    %     obj.interfaceList = sort([obj.interfaceList interfId]);
    %     obj.interfaces{end+1} = interf;
    %   end
    %end


    function assembleSystem(obj,dt)
      % loop trough solver database and compute non-costant jacobian
      % blocks and rhs block
      for solver = obj.solverNames
        assembleSystem(obj.getPhysicsSolver(solver),dt);
      end
    end


    function addPhysicsSolver(obj,solverInput)

      assert(isempty(getJacobian(obj)),"Cannot add a physics solver " + ...
        "after system has already been assembled");

      % Add a new solver to the Discretizer
      assert(nargin == 2,"Input must be an xml file or a scalar struct")

      if ~isstruct(solverInput)
        solverInput = readstruct(solverInput,AttributeSuffix="");
      end

      if isfield(solverInput,"Solver")
        solverInput = solverInput.Solver;
      end

      obj.solverNames = string(fieldnames(solverInput));
      obj.solverNames = reshape(obj.solverNames,1,[]);

      for solverName = obj.solverNames
        % create and register the solver
        solver = feval(solverName,obj);
        solver.registerSolver(solverInput.(solverName));
        obj.physicsSolvers(solverName) = solver;
      end

      nV = obj.dofm.getNumberOfVariables();
      obj.J = cell(nV);
      obj.rhs = cell(nV,1);

    end



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
        J = obj.J(id,id);
      elseif nargin == 3
        idRow = obj.dofm.getVariableId(varargin{1});
        idCol = obj.dofm.getVariableId(varargin{2});
        J = obj.J(idRow,idCol);
      else
        error("Too many input arguments")
      end

    end


    function [Jum,Jmu] = getInterfaceJacobian(obj,interfaceId)

      % get the multidomain coupling blocks for a certain interface
      Jum = obj.Jum{interfaceId};
      Jmu = obj.Jmu{interfaceId};

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
        rhs = obj.rhs(id);
      else
        error("Too many input arguments")
      end
    end

    % function printState(obj)
    %   % print solution of the model according to the print time in the
    %   % list
    %
    %   % if obj.state.t >= obj.simparams.tMax
    %   %
    %   %   time = getState(obj).t;
    %   %   writeVTK(obj,time);
    %   %   writeMatFile(obj,time,obj.outstate.timeID);
    %   %
    %   % else
    %
    %   if obj.outstate.timeID <= length(obj.outstate.timeList)
    %
    %     time = obj.outstate.timeList(obj.outstate.timeID);
    %
    %     % loop over print times within last time step
    %     while time <= obj.state.t
    %
    %       assert(time >= obj.stateOld.t, 'Print time %f out of range (%f - %f)',...
    %         time, obj.stateOld.t, obj.state.t);
    %
    %       assert(obj.state.t - obj.stateOld.t > eps('double'),...
    %         'Time step is too small for printing purposes');
    %
    %       % move this into the solution scheme class
    %       % compute factor to interpolate current and old state variables
    %       fac = (time - obj.stateOld.t)/(obj.state.t - obj.stateOld.t);
    %       if isnan(fac) || isinf(fac)
    %         fac = 1;
    %       end
    %
    %       writeVTK(obj,fac,time);
    %
    %       writeMatFile(obj,fac,obj.outstate.timeID);
    %
    %       obj.outstate.timeID = obj.outstate.timeID + 1;
    %
    %       if obj.outstate.timeID > length(obj.outstate.timeList)
    %         break
    %       else
    %         time = obj.outstate.timeList(obj.outstate.timeID);
    %       end
    %
    %     end
    %
    %   end
    % end


    function out = writeVTK(obj,fac,time)
      % write results to VTKoutput

      obj.vtmBlock = obj.outstate.vtkFile.createElement('Block');

      cellData3D = struct('name', [], 'data', []);
      pointData3D = struct('name', [], 'data', []);

      for solv = obj.solverNames
        solver = getPhysicsSolver(obj,solv);
        [cellData,pointData] = writeVTK(solver,fac,time);
        cellData3D = OutState.mergeOutFields(cellData3D,cellData);
        pointData3D = OutState.mergeOutFields(pointData3D,pointData);
      end

      cellData3D = OutState.printMeshData(obj.grid.topology,cellData3D);

      % write dataset to vtmBlock
      obj.outstate.writeVTKfile(obj.vtmBlock,obj.getOutName(),obj.grid.topology,...
        time, pointData3D, cellData3D, [], []);

      out = obj.vtmBlock;

    end


    function writeMatFile(obj,fac,timeID)
      % write to MAT-file
      obj.outstate.matFile(timeID).time = obj.outstate.timeList(timeID);

      for solv = obj.solverNames
        getPhysicsSolver(obj,solv).writeMatFile(fac,timeID);
      end
    end


    function nDofs = getNumbDoF(obj,varargin)

      nDofs = obj.dofm.getNumbDoF(varargin{:});

    end
    %end

  end

  methods(Access = private)

    function setDiscretizer(obj,varargin)

      % key-value input
      setInput(obj,varargin{:});

      % validate consistency between input
      validateInput(obj);

      obj.state = State();

      if ~isempty(obj.grid.topology)
        obj.dofm = DoFManager(obj.grid.topology);
      end

    end

    function setInput(obj, varargin)

      msg = 'Invalid key-value pair for Discretizer class \n';

      % Check that we have an even number of inputs
      if mod(length(varargin), 2) ~= 0
        error('Arguments must come in key-value pairs.');
      end

      % Loop through the key-value pairs
      for k = 1:2:length(varargin)
        key = varargin{k};
        value = varargin{k+1};

        if isempty(value)
          continue
        end

        if ~ischar(key) && ~isstring(key)
          error('Keys must be strings');
        end

        switch lower(key)
          % case 'simulationparameters'
          %   assert(isa(value, 'SimulationParameters')|| isempty(value),msg)
          %   obj.simparams = value;
          case 'grid'
            obj.grid = value;
            % check that grid has been defined correctly
            if ~isfield(obj.grid,"topology")
              obj.grid.topology = [];
            end
            if ~isfield(obj.grid,"cells")
              obj.grid.cells = [];
            end
            if ~isfield(obj.grid,"faces")
              obj.grid.faces = [];
            end

            isGridCorrect = numel(fieldnames(obj.grid))==3;

            assert(isGridCorrect,"Error in Discretizer: grid input is not correct. " + ...
              "grid must be a struct with fields: 'topology','cells','faces'.");

          case 'materials'
            assert(isa(value, 'Materials'),msg)
            obj.materials = value;
          case 'boundaries'
            assert(isa(value, 'Boundaries'),msg)
            obj.bcs = value;
          otherwise
            error('Unknown input key %s for Discretier \n', key);
        end
      end


    end

    function validateInput(obj)
      msh = obj.grid.topology;

      u = unique(msh.cellTag);
      if ~isempty(u)
        assert(max(u)==length(u),"cellTag numbering must be progressive from 1 to the number of cellTags")
        msh.nCellTag = max(u);
      end

      u = unique(msh.surfaceTag);
      if ~isempty(u)
        assert(max(u)==length(u),"surfaceTag numbering must be progressive from 1 to the number of cellTags")
        msh.nSurfaceTag = max(u);
      end

    end


    function outName = getOutName(obj)

      outName = sprintf('Domain_%i',obj.domainId);


    end

  end

end

