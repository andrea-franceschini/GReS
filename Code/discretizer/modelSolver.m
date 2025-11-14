classdef modelSolver < handle
  % General model class
  % This class stores all the information related to the problem at hand
  % and store a map to the existing domainSolver and interfaceSolver



  properties (GetAccess=public, SetAccess=private)
    % model commons
    physicsSolvers         % db for physics solver in the model
    interfaceSolvers      % db for interface solvers in the model
    simparams             % simulation parameters
  end


  methods (Access = public)
    function obj = modelSolver(simparams)

      obj.physicsSolvers = containers.Map('KeyType','double','ValueType','any');
      obj.interfaceSolvers = containers.Map('KeyType','double','ValueType','any');

    end


    function runProblem()
    end

    function printSolution()
    end


    function out = getSolver(obj,id)
      % varargin: string 1 or 2 field identifier
      out = obj.physicsSolvers(id);
    end

    function out = getInterface(obj,id)
      % varargin: string 1 or 2 field identifier
      out = obj.interfaceSolvers(id);
    end


    function addInterface(obj,interfId,interf)
      % add mortar interface to current domain
      if ~ismember(interfId,obj.interfaceList)
        obj.interfaceList = sort([obj.interfaceList interfId]);
        obj.interfaces{end+1} = interf;
      end
    end


    function computeMatricesAndRhs(obj,stateOld,dt)
      % loop trough solver database and compute non-costant jacobian
      % blocks and rhs block
      for i = 1:numel(obj.solverNames)
        computeMat(obj.solver(obj.solverNames{i}),stateOld,dt);
        computeRhs(obj.solver(obj.solverNames{i}),stateOld,dt);
      end
    end

    function out = isCoupled(obj,field1,field2)
      % check if input fields are coupled, i.e exist cells having both
      % fields activated
      sub1 = getActiveSubdomain(obj.dofm,field1);
      sub2 = getActiveSubdomain(obj.dofm,field2);
      out = any(intersect(sub1,sub2));
    end

    function initState(obj)
      % loop trough active single physics solver and update the state class
      % accordingly
      for i = 1:numel(obj.fields)
        % loop trough active fields and update the state structure
        initState(obj.getSolver(obj.fields(i)));
      end
    end

    function updateState(obj,du)
      % update current state
      for i = 1:numel(obj.fields)
        obj.getSolver(obj.fields(i)).updateState(du);
      end
    end
  end

  methods(Access = private)
    function setDiscretizer(obj,varargin)

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
          fieldNames{k} = Discretizer.getFieldString(flds(i),flds(j));
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
          case 'modeltype'
            assert(isa(value, 'ModelType')|| isempty(value),msg)
            obj.model = value;
          case 'simulationparameters'
            assert(isa(value, 'SimulationParameters')|| isempty(value),msg)
            obj.simparams = value;
          case 'dofmanager'
            assert(isa(value, 'DoFManager') || isempty(value),msg)
            obj.dofm = value;
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

    function checkTimeDependence(obj)
      if isempty(obj.simparams)
        return
      end
      % check if there is any time dependence in the input model
      % this must be moved into the single physics models
      if ~isSinglePhaseFlow(obj.model)
        % Biot model is time dependent
        setTimeDependence(obj.simparams,false);
        return
      else
        % check if fluid is incompressible
        beta = getFluidCompressibility(obj.materials.getFluid());
        if beta < eps
          setTimeDependence(obj.simparams,false);
        end
      end
    end

    function addPhysics(obj,f1,f2,state)
      % Add new key to solver database
      % Prepare input fields for solver definition
      if ~isCoupled(obj,f1,f2)
        return
      end

      fList = Discretizer.getFieldString(f1,f2);

      if ~obj.solverMap.isKey(fList)
        error('A physical module coupling %s with %s is not available.',f1,f2)
      else
        solv = obj.solverMap(fList);
      end

      obj.solver(fList) = solv(obj.model,obj.simparams,obj.dofm,...
        obj.grid,obj.materials,obj.bcs,state);
    end


    function setSolverMap(obj)

      obj.solverMap = containers.Map('KeyType','char','ValueType','any');

      subClasses = [findSubClasses('SinglePhysics','SinglePhysics'), ...
        findSubClasses('CouplingPhysics','CouplingPhysics')];

      for i = 1:numel(subClasses)
        obj.solverMap = feval([subClasses{i} '.registerSolver'],...
          obj.solverMap,subClasses{i});
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

    function f = getFieldString(varargin)
      if numel(varargin) == 1
        f = char(varargin{1});
      elseif numel(varargin) == 2
        f = join(unique(sort({char(varargin{1}),char(varargin{2})})));
        f = f{:};
      else
        error('Too many input fields')
      end
    end

    function mat_new = expandMat(mat,n)
      % Get the size of the original matrix
      [s1, s2] = size(mat);

      % Initialize the sparse matrix: row indices, column indices, and values
      rows = [];
      cols = [];
      values = [];

      % Loop through the original matrix and populate the sparse matrix
      for s = n-1:-1:0
        % Get the row and column indices for the block
        r1 = n*(1:s1) - s;
        r2 = n*(1:s2) - s;
        [colIdx,rowIdx]  = meshgrid(r2,r1);
        % Add the values from the original matrix to the sparse matrix
        rows = [rows; rowIdx(:)];
        cols = [cols; colIdx(:)];
        values = [values; mat(:)];
      end

      % Create the sparse matrix directly from the row, column indices, and values
      mat_new = sparse(rows, cols, values, s1 * n, s2 * n);
    end
  end
end
