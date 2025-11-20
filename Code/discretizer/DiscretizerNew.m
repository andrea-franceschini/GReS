classdef DiscretizerNew < handle
  % General discretizer class
  properties (GetAccess=public, SetAccess=private)

    physicsSolvers      % physics solvers database
    dofm                % dofManager
    simparams
    bcs
    outstate
    materials
    solverNames
    grid = struct('topology',[],'faces',[],'cells',[])

  end

  properties (GetAccess=private, SetAccess=public)
    interfaceList = [];
    interfaces = []
    interfaceSurf            % surfaceTag in each surface of interfaceSurf
    % empty - single domain simulation
    % not empty - call to mesh glue instances

    state
    stateOld
  end

  properties (GetAccess=private, SetAccess=public)
    J
    rhs
  end


  methods (Access = public)
    function obj = DiscretizerNew(varargin)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      obj.physicsSolvers = containers.Map('KeyType','char','ValueType','any');
      obj.setDiscretizer(varargin{:});
    end

    function applyBC(obj,t)
      bcList = keys(obj.bcs.keys);

      for bcId = string(bcList)
        % loop over available bcs
        bcVar = obj.bcs.getVariable(bcId);

        for solv = obj.solverNames
          % loop over available solvers
          solver = obj.getPhysicsSolver(solv);

          % check if BCs applies to the solver
          if any(strcmpi(solver.getField(),bcVar))
            solver.applyBC(t,bcId,bcVar)
          end

        end
      end

    end


    function applyDirVal(obj,t)
      bcList = keys(obj.bcs.keys);

      for bcId = string(bcList)
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
            solver.applyBDirVal(t,bcId,bcVar)
          end
        end
      end
    end


    function out = getPhysicsSolver(obj,id)
      out = obj.physicsSolvers(id);
    end

    function out = getState(obj)
      out = obj.state;
    end

    function out = getStateOld(obj)
      out = obj.stateOld;
    end


    function printState(obj,t)
      % print solution of the model according to the print time in the
      % list

      if t >= obj.simparams.tMax

        time = obj.state.t;

        writeVTK(obj,time);

        writeMatFile(obj,time);

        % loop over print times
      else

        if obj.outstate.timeID <= length(obj.outstate.timeList)

          time = obj.outstate.timeList(obj.outstate.timeID);

          while time <= obj.state.t

            assert(time > obj.stateOld.t, 'Print time %f out of range (%f - %f)',...
              time, obj.stateOld.t, stateNew.t);

            assert(obj.state.t - obj.stateOld.t > eps('double'),...
              'Time step is too small for printing purposes');

            writeVTK(obj,time);

            writeMatFile(obj,time);

            obj.outstate.timeID = obj.outstate.timeID + 1;
            
            time = obj.outstate.timeList(obj.outstate.timeID);

          end

        end
      end
    end
 





    function out = getSolver(obj,varargin)
      % varargin: string 1 or 2 field identifier
      fld = DiscretizerNew.getFieldString(varargin{:});
      out = obj.solver(fld);
    end


    function addInterface(obj,interfId,interf)
      % add mortar interface to current domain
      if ~ismember(interfId,obj.interfaceList)
        obj.interfaceList = sort([obj.interfaceList interfId]);
        obj.interfaces{end+1} = interf;
      end
    end


    function assembleSystem(obj,dt)
      % loop trough solver database and compute non-costant jacobian
      % blocks and rhs block
      for i = 1:numel(obj.solverNames)
        assembleSystem(obj.solver(obj.solverNames{i}),dt);
      end
    end


  end

  methods(Access = private)
    function setDiscretizer(obj,varargin)

      % key-value input 
      setInput(obj,varargin{:});


      obj.setSolverMap();
      flds = getFieldList(obj.dofm);
      nF = numel(flds);
      % loop over all fields and define corresponding models
      k = 0;
      % create the handle to state object that will be shared across all physical
      % modules
      stat = State();
      fieldNames = cell(nF^2,1);
      for i = 1:nF
        for j = i:nF
          k = k+1;
          addPhysics(obj,flds(i),flds(j),stat);
          fieldNames{k} = DiscretizerNew.getFieldString(flds(i),flds(j));
        end
      end
      obj.state = stat;
      obj.fields = flds;
      obj.solverNames = fieldNames(1:k);
      obj.interfaces = cell(0);
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

        if ~ischar(key) && ~isstring(key)
          error('Keys must be strings');
        end

        switch lower(key)
          case 'simulationparameters'
            assert(isa(value, 'SimulationParameters')|| isempty(value),msg)
            obj.simparams = value;
          case 'grid'
            assert(isstruct(value),msg)
            obj.grid = value;
          case 'materials'
            assert(isa(value, 'Materials') || isempty(value),msg)
            obj.materials = value;
          case 'boundaries'
            assert(isa(value, 'Boundaries') || isempty(value),msg)
            obj.bcs = value;
          case 'outstate'
            assert(isa(value, 'OutState') || isempty(value),msg)
            obj.outstate = value;
          otherwise
            error('Unknown input %s for Discretier \n', key);
        end
      end

      % check that grid has been defined
      assert(~isempty(obj.grid),['Grid structure with a topology field ' ...
        'is required for Discretizer class']);
    end


    function addPhysicsSolver(obj,solverInput)

      % Add a new solver to the Discretizer
      assert(nargin == 2,"Input must be an xml file or a scalar struct")

      if isstruct(solverInput)
        solverInput = readstruct(solverInput,AttributeSuffix="");
      end

      if isfield(solverInput,"Solver")
        solverInput = solverInput.Solver;
      end

      obj.solverNames = string(fieldnames(solverInput));

      for i = 1:numel(obj.solverNames)
        % create and register the solver
        solvName = obj.solverNames(i);
        solv = feval(solvName,obj);
        solv.registerSolver(solverInput.(solvName))
        obj.physicsSolvers(solvName)
      end


    end


    function writeVTK(obj,time)

      if obj.outstate.writeVtk

        for solv = obj.solverNames
          solver = getPhysicsSolver(obj,solv);

          cellData3D = [];
          pointData3D = [];
          [cellData,pointData] = printVTK(solver,time);
          cellData3D = OutState.mergeOutFields(cellData3D,cellData);
          pointData3D = OutState.mergeOutFields(pointData3D,pointData);
        end

        cellData3D = OutState.printMeshData(obj.grid.topology,cellData3D);
        obj.outstate.VTK.writeVTKFile(time, pointData3D, cellData3D, [], []);
      end

    end

    function writeMatFile(obj,time)
      
      % write to MAT-file
      if obj.outstate.writeSolution

        tID = obj.outstate.timeID;

        for solv = obj.solverNames
          solver = getPhysicsSolver(obj,solv);
          writeMatFile(solver,tID,time);
        end

      end
    end


  end

  methods (Static)
    %      function [row,col,val,c] = computeLocalMatrix(mat,row,col,val,c,w,dofRow,dofCol)
    %        % shortcut for assemblying local matrix contributions in sparse format
    %        mat = mat.*reshape(w,1,1,[]);
    %        mat = sum(mat,3);
    %        n = numel(mat);
    %        [J, I] = meshgrid(1:size(mat,2), 1:size(mat,1));
    %        row(c+1:c+n) = dofRow(I);
    %        col(c+1:c+n) = dofCol(J);
    %        val(c+1:c+n) = mat(:);
    %        c = c+n;
    %      end

    % 
    % 
    % function mat_new = expandMat(mat,n)
    %   % Get the size of the original matrix
    %   [s1, s2] = size(mat);
    % 
    %   % Initialize the sparse matrix: row indices, column indices, and values
    %   rows = [];
    %   cols = [];
    %   values = [];
    % 
    %   % Loop through the original matrix and populate the sparse matrix
    %   for s = n-1:-1:0
    %     % Get the row and column indices for the block
    %     r1 = n*(1:s1) - s;
    %     r2 = n*(1:s2) - s;
    %     [colIdx,rowIdx]  = meshgrid(r2,r1);
    %     % Add the values from the original matrix to the sparse matrix
    %     rows = [rows; rowIdx(:)];
    %     cols = [cols; colIdx(:)];
    %     values = [values; mat(:)];
    %   end
    % 
    %   % Create the sparse matrix directly from the row, column indices, and values
    %   mat_new = sparse(rows, cols, values, s1 * n, s2 * n);
    % end
  end
end
