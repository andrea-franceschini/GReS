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
      % obj.bcs = Boundaries();
      % obj.outstate = OutState();
      % obj.materials = Materials();
      % obj.grid = Grid();

      obj.physicsSolvers = containers.Map('KeyType','char','ValueType','any');
      obj.setDiscretizer(varargin{:});
    end


    function applyBC(obj,t,varargin)
      bcList = obj.bcs.getBCList();

      if isempty(varargin)
        varNames = obj.dofm.getVariableNames();
      else
        varNames = [varargin{:}];
      end

      for i = 1:numel(bcList)
        bcId = bcList(i);

        % loop over available bcs
        bcVar = obj.bcs.getVariable(bcId);

        if ~strcmp(bcVar,varNames)
          continue
        end

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

      for i = 1:numel(bcList)
        bcId = bcList(i);
        
        % discard non-dirichlet BC
        if ~isEssential(obj.bcs,bcId)
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


    function addPhysicsSolvers(obj,solverInput)

      assert(isempty(getJacobian(obj)),"Cannot add a physics solver " + ...
        "after system has already been assembled");

      % Add a new solver to the Discretizer
      assert(nargin == 2,"Input must be an xml file or a struct")

      solverInput = readInput(solverInput);
      if isfield(solverInput,"Solver")
        solverInput = solverInput.Solver;
      end

      sN = string(fieldnames(solverInput));
      sN = reshape(sN,1,[]);

      % assert(numel(sN)==numel(solverInput),"A PhysicsSolver cannot" + ...
      %   "be defined more than once in a model.")

      for solverName = sN
        % create and register the solver
        addPhysicsSolver(obj,solverName,solverInput.(solverName));
      end

      % nV = obj.dofm.getNumberOfVariables();
      % obj.J = cell(nV);
      % obj.rhs = cell(nV,1);

    end

    function addPhysicsSolver(obj,solverName,varargin)

      if ~any([ischar(solverName),isstring(solverName)])
        error("First input must be the name of the PhysicsSolver")
      end

      if any(strcmp(solverName,obj.solverNames))
        error('Solver %s has already been defined',solverName)
      else
        obj.solverNames = [obj.solverNames, string(solverName)];
      end
      solver = feval(solverName,obj);
      solver.registerSolver(varargin{:});
      obj.physicsSolvers(solverName) = solver;

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


    function initialize(obj)

      % prepare the discretizer before starting a new simulation

      % initialize block jacobian and rhs
      nV = obj.dofm.getNumberOfVariables();
      obj.J = cell(nV);
      obj.rhs = cell(nV,1);

      prepareBoundaryConditions(obj);

      for solver = obj.solverNames
        initialize(obj.getPhysicsSolver(solver));
      end

    end

    function timeStepSetup(obj)

      % prepare the physics for each new time step
      for solver = obj.solverNames
        timeStepSetup(obj.getPhysicsSolver(solver));
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

      cellData3D = OutState.printMeshData(obj.grid,cellData3D);

      % write dataset to vtmBlock
      obj.outstate.writeVTKfile(obj.vtmBlock,obj.getOutName(),obj.grid,...
        time, pointData3D, cellData3D, [], []);

      out = obj.vtmBlock;

    end


    function writeSolution(obj,fac,timeID)
      % write to MAT-file
      obj.outstate.results(timeID).time = obj.outstate.timeList(timeID);

      for solv = obj.solverNames
        getPhysicsSolver(obj,solv).writeSolution(fac,timeID);
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

      if ~isempty(obj.grid)
        obj.dofm = DoFManager(obj.grid);
      end

    end

    function setInput(obj, varargin)

      default = struct('grid',Grid(),...
                       'materials',Materials(),...
                       'boundaries',Boundaries());

      % make read case-insensitive
      input = readInput(varargin{:});
      fn = fieldnames(input);
      for i = 1:numel(fn)
        newName = lower(fn{i});
        if ~strcmp(fn{i}, newName)
          input.(newName) = input.(fn{i});
          input = rmfield(input, fn{i});
        end
      end

      params = readInput(default,input);

      obj.grid = params.grid;
      obj.materials = params.materials;
      obj.bcs = params.boundaries;

    end

   
    function validateInput(obj)

      if obj.grid.nNodes == 0
        % the grid is still empty, nothing to validate
        return
      end

      msh = obj.grid.cells;

      t = unique(msh.tag);
      if ~isempty(t)
        assert(t(end) == length(t),"cellTag numbering must be progressive from 1 to the number of cellTags")
      end
      msh.nTag = t(end);
      obj.grid.cells = msh;

      msh = obj.grid.surfaces;

      t = unique(msh.tag);
      if ~isempty(t)
        assert(t(end) == length(t),"surfaceTag numbering must be progressive from 1 to the number of cellTags")
      end
      msh.nTag = t(end);
      obj.grid.surfaces = msh;

    end

    function prepareBoundaryConditions(obj)

      % preprocess the boundary condition once the type of the target field
      % is knwon

      setBCList(obj.bcs);

      bcList = obj.bcs.getBCList();

      for i = 1:numel(bcList)
        bcId = bcList(i);
        
        % loop over available bcs
        bcVar = obj.bcs.getVariable(bcId);

        targetField = obj.dofm.getFieldLocation(bcVar);

        obj.bcs.initialize(bcId,targetField);

      end
 
    end



    function outName = getOutName(obj)
       % property domainId not set yet
      outName = sprintf('Domain_%i',obj.domainId);


    end

  end


  methods (Static)


    function checkClass(val,fldName,className)

       if ~isa(val, className)
         error("Invalid key-value pair for Discretizer class. Input field %s must be an object of class %s \n",fldName,className);
       end
    end


  end

end

