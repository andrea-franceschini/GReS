classdef Discretizer < handle
  % General discretizer class
  properties (GetAccess=public, SetAccess=private)
    solver      % database for physics solvers in the model
    dofm        % dofManager
    simparams
    bcs
    outstate
    materials
    solverNames
    model
    fields
    grid
  end

  properties (GetAccess=public, SetAccess=public)
    interfaceList = [];
    interfaces = []
    interfaceSurf            % surfaceTag in each surface of interfaceSurf
    % empty - single domain simulation
    % not empty - call to mesh glue instances
    state
  end

  properties (Access = private)
    solverMap
  end

  methods (Access = public)
    function obj = Discretizer(varargin)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      obj.solver = containers.Map('KeyType','char','ValueType','any');
      obj.setDiscretizer(varargin{:});
      obj.initState();
      obj.checkTimeDependence();
    end

    function applyBC(obj,t)
      % Apply boundary condition to blocks of physical solver
      % ents: id of constrained entity
      % vals: value to assign to each entity
      bcList = obj.bcs.db.keys;
      % get entities and values of boundary condition
      for bcId = string(bcList)
        field = obj.bcs.getPhysics(bcId);
        % get id of constrained entities and corresponding BC values
        [bcEnts,bcVals] = getBC(getSolver(obj,field),bcId,t);

        %removeInterfaceBC(obj,bcEnts,bcVals)
        % apply Boundary conditions to each Jacobian/rhs block
        for f = obj.fields
          if ~isCoupled(obj,field,f)
            continue
            % skip pair of uncoupled physics
          end
          switch obj.bcs.getType(bcId)
            case {'Dirichlet','Seepage'}
              % if nargin > 2
              %   assert(~isempty(obj.interfaceList),['Too many input arguments: ' ...
              %     'invalid domain id input for single domain BC imposition']);
              %   for i = 1:length(obj.interfaceList)
              %     [bcEnts,bcVals] = obj.interfaces{i}.removeSlaveBCdofs(field,[bcEnts,bcVals],idDom);
              %   end
              % end
              applyDirBC(obj.getSolver(field,f),field,bcEnts,bcVals);
            case {'Neumann','VolumeForce'}
              applyNeuBC(obj.getSolver(field,f),bcEnts,bcVals);
          end
        end
      end
    end

    function applyDirVal(obj,t)
      % Apply boundary condition to blocks of physical solver
      % ents: id of constrained entity
      % vals: value to assign to each entity
      bcList = obj.bcs.db.keys;
      % get entities and values of boundary condition
      for bcId = string(bcList)
        if ~strcmp(obj.bcs.getType(bcId),'Dirichlet')
          continue
        end
        field = obj.bcs.getPhysics(bcId);
        [bcEnts,bcVals] = getBC(getSolver(obj,field),bcId,t);
        % if nargin > 2
        %   assert(~isempty(obj.interfaceList),['Too many input arguments: ' ...
        %     'invalid domain id input for single domain BC imposition']);
        %   for i = 1:length(obj.interfaceList)
        %     [bcEnts,bcVals] = obj.interfaces{i}.removeSlaveBCdofs(field,[bcEnts,bcVals],idDom);
        %   end
        % end
        getSolver(obj,field).applyDirVal(bcEnts,bcVals);
      end
    end

    function J = assembleJacobian(obj)
      % put together jacobian blocks of SinglePhysicsSolver and
      % CoupledSolver in the model
      % Use the ordering specified in DoF manager class
      switch obj.dofm.ordering
        case 'field'
          nFld = numel(obj.fields);
          J = cell(nFld,nFld);
          for i = 1:nFld
            for j = 1:nFld
              J{i,j} = getJacobian(obj.getSolver(obj.fields(i),obj.fields(j)),obj.fields(i));
            end
          end
        otherwise
          error('Invalid DoF manager ordering')
      end
    end

    function rhs = assembleRhs(obj)
      % put together rhs blocks of SinglePhysicsSolver and
      % CoupledSolver in the model
      nFld = numel(obj.fields);
      rhs = cell(nFld,1);
      for i = 1:nFld
        rhs{i} = zeros(getNumDoF(obj.dofm,obj.fields(i)),1);
        for j = 1:nFld
          rhs{i} = rhs{i} + ...
            getRhs(getSolver(obj,obj.fields(i),obj.fields(j)),obj.fields(i));
        end
      end
    end

     function printState(obj,stateOld)
      % print solution of the model according to the print time in the
      % list
      % Initialize input structure for VTK output
      cellData3D = [];
      pointData3D = [];
      if nargin == 1
        time = obj.state.t;
        % print result to mat-file
        % this is not modular and will be updated in future version of the code
        if obj.outstate.writeSolution
          % obj.outstate.results.expTime(obj.outstate.timeID,1) = time;
          obj.outstate.results(obj.outstate.timeID).expTime = time;
          if isPoromechanics(obj.model)
            % obj.outstate.results.expDispl(:,obj.outstate.timeID) = obj.state.data.dispConv;
            obj.outstate.results(obj.outstate.timeID).expDispl = obj.state.data.dispConv;
          end
          if isFlow(obj.model)
            % obj.outstate.results.expPress(:,obj.outstate.timeID) = obj.state.data.pressure;
            obj.outstate.results(obj.outstate.timeID).expPress = obj.state.data.pressure;
          end
          if isVariabSatFlow(obj.model)
            % obj.outstate.results.expSat(:,obj.outstate.timeID) = obj.state.data.saturation;
            obj.outstate.results(obj.outstate.timeID).expSat = obj.state.data.saturation;
          end
        end
        for fld = obj.fields
          [cellData,pointData] = printState(obj.getSolver(fld),obj.state);
          cellData3D = [cellData3D; cellData];
          pointData3D = [pointData3D; pointData];
        end
        cellData3D = OutState.printMeshData(obj.grid.topology,cellData3D);
        if obj.outstate
          obj.outstate.VTK.writeVTKFile(time, pointData3D, cellData3D, [], []);
        end
        % update the print structure
      elseif nargin == 2
        stateNew = obj.state;
        if obj.outstate.timeID <= length(obj.outstate.timeList)
          while (obj.outstate.timeList(obj.outstate.timeID) <= stateNew.t)
            assert(obj.outstate.timeList(obj.outstate.timeID) > stateOld.t, ...
              'Print time %f out of range (%f - %f)',obj.outstate.timeList(obj.outstate.timeID), ...
              stateOld.t,stateNew.t);
            assert(stateNew.t - stateOld.t > eps('double'),'Dt too small for printing purposes');
            %
            time = obj.outstate.timeList(obj.outstate.timeID);
            if obj.outstate.writeSolution
              % print solution to mat-file
              fac = (time - stateOld.t)/(stateNew.t - stateOld.t);
              % obj.outstate.results.expTime(obj.outstate.timeID+1,1) = time;
              obj.outstate.results(obj.outstate.timeID+1).expTime = time;
              if isPoromechanics(obj.model)
                % obj.outstate.results.expDispl(:,obj.outstate.timeID+1) = stateNew.data.dispConv*fac+stateOld.data.dispConv*(1-fac);
                obj.outstate.results(obj.outstate.timeID+1).expDispl = stateNew.data.dispConv*fac+stateOld.data.dispConv*(1-fac);
              end
              if isFlow(obj.model)
                % obj.outstate.results.expPress(:,obj.outstate.timeID+1) = stateNew.data.pressure*fac+stateOld.data.pressure*(1-fac);
                obj.outstate.results(obj.outstate.timeID+1).expPress = stateNew.data.pressure*fac+stateOld.data.pressure*(1-fac);
              end
              if isVariabSatFlow(obj.model)
                % obj.outstate.results.expSat(:,obj.outstate.timeID+1) = stateNew.data.saturation*fac+stateOld.data.saturation*(1-fac);
                obj.outstate.results(obj.outstate.timeID+1).expSat = stateNew.data.saturation*fac+stateOld.data.saturation*(1-fac);
              end
            end
            % Write output structure looping through available models
            for fld = obj.fields
              [cellData,pointData] = printState(obj.getSolver(fld),...
                stateOld, stateNew, time);
              % merge new fields
              cellData3D = OutState.mergeOutFields(cellData3D,cellData);
              pointData3D = OutState.mergeOutFields(pointData3D,pointData);
            end
            cellData3D = OutState.printMeshData(obj.grid.topology,cellData3D);
            if obj.outstate.writeVtk
              obj.outstate.VTK.writeVTKFile(time, pointData3D, cellData3D, [], []);
            end
            obj.outstate.timeID = obj.outstate.timeID + 1;
            if obj.outstate.timeID > length(obj.outstate.timeList)
              break
            end
          end
        end
      end
    end


    function out = getSolver(obj,varargin)
      % varargin: string 1 or 2 field identifier
      fld = Discretizer.getFieldString(varargin{:});
      out = obj.solver(fld);
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
