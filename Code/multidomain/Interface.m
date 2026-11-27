classdef Interface < handle
  % Geometry, state, and system blocks shared by interface constraints.

  properties
    model
    grid
    outstate
    state
    id
    dofm
    vtmBlock
  end

  properties (SetAccess=private)
    domains
    quadrature
    constraints
    constraintNames = strings(1,0)
  end

  properties (Access=private)
    initialized = false
  end

  methods
    function obj = Interface(input,varargin)
      obj.domains = getInterfaceGrid(obj,input,varargin{:});
      obj.dofm = DoFManager(obj.grid);
      obj.state = State();
      obj.constraints = containers.Map('KeyType','char','ValueType','any');

    end

    function setMortarQuadrature(obj,multType,varargin)

      if ~isempty(obj.quadrature)
        return
      end

      % every interface can have one single quadrature object
      % it is not possible to use different quadrature rules on different
      % constraints

      assert(obj.dofm.getNumberOfVariables() == 0,...
        'Set up mortar quadrature before registering interface fields');
      options = readInput(struct('type',"SegmentBasedQuadrature"),varargin{:});
      quad = feval(options.type,obj.grid,multType,options);
      quad.processMortarTriplets();

      % Pruning can replace meshes with subgrids that have no edge topology.
      grids = {obj.grid,obj.grid.traceGrids(1),obj.grid.traceGrids(2)};
      for i = 1:numel(grids)
        mesh = grids{i};
        if mesh.nDim == 2 && mesh.faces.num == 0
          mesh.processFaces();
        end
      end

      obj.quadrature = quad;

    end


    function addConstraints(obj,input)
      input = readInput(input);
      if isfield(input,"Constraint")
        input = input.Constraint;
      end

      for name = reshape(string(fieldnames(input)),1,[])
        obj.addConstraint(name,input.(name));
      end
    end

    function addConstraint(obj,name,varargin)
      
      assert(~obj.initialized, ...
        "Cannot add a constraint after the interface has been initialized");
      assert((ischar(name) && isrow(name)) || ...
        (isstring(name) && isscalar(name) && ~ismissing(name)), ...
        "First input must be the name of the Constraint");
      name = char(name);
      assert(~isempty(name),"Constraint name must be nonempty");
      assert(~isKey(obj.constraints,name), ...
        "Constraint %s has already been defined",name);

      constraint = feval(name,obj);
      assert(isa(constraint,'Constraint'),"%s must derive from Constraint",name);
      constraint.registerConstraint(varargin{:});
      obj.constraints(name) = constraint;
      obj.constraintNames(end+1) = string(name);
    end

    function constraint = getConstraint(obj,name)
      constraint = obj.constraints(char(name));
    end

    function names = getConstraintNames(obj)

      names = obj.constraintNames;
    end

    function initialize(obj)
      for name = obj.constraintNames
        obj.getConstraint(name).initialize();
      end
      obj.initialized = true;
    end

    function timeStepSetup(obj)
      for name = obj.constraintNames
        obj.getConstraint(name).timeStepSetup();
      end
    end

    function assembleSystem(obj,varargin)
      for name = obj.constraintNames
        obj.getConstraint(name).assembleConstraint(varargin{:});
      end
    end

    function updateState(obj,solution)
      for name = obj.constraintNames
        obj.getConstraint(name).updateState(solution);
      end
    end

    function hasConfigurationChanged = updateConfiguration(obj)
      hasConfigurationChanged = false;
      for name = obj.constraintNames
        changed = obj.getConstraint(name).updateConfiguration();
        hasConfigurationChanged = hasConfigurationChanged || changed;
      end
    end

    function resetConfiguration(obj)
      for name = obj.constraintNames
        obj.getConstraint(name).resetConfiguration();
      end
    end

    function advanceState(obj)
      for name = obj.constraintNames
        obj.getConstraint(name).advanceState();
      end
    end

    function goBackState(obj)
      for name = obj.constraintNames
        obj.getConstraint(name).goBackState();
      end
    end

    function out = isLinear(obj)
      out = true;
      for name = obj.constraintNames
        out = out && obj.getConstraint(name).isLinear();
      end
    end

    function nDofs = getNumbDoF(obj,varargin)
      nDofs = obj.dofm.getNumbDoF(varargin{:});
    end

    function dofm = getDoFManager(obj,side)
      if nargin == 1
        dofm = obj.dofm;
      else
        dofm = obj.domains(side).dofm;
      end
    end

    function val = getVariableField(obj,side,domainField)
      domain = obj.domains(side);
      ents = domain.dofm.getActiveEntities(domainField,1);
      val = domain.getState(domainField);
      val = val(ents);
    end

    function block = writeVTK(obj,fac,time)
      obj.vtmBlock = obj.outstate.vtkFile.createElement('Block');
      cellData = struct('name', [], 'data', []);
      pointData = struct('name', [], 'data', []);

      for name = obj.constraintNames
        [constraintCells,constraintPoints] = obj.getConstraint(name).writeVTK(fac,time);
        cellData = OutState.mergeOutFields(cellData,constraintCells);
        pointData = OutState.mergeOutFields(pointData,constraintPoints);
      end

      cellData = OutState.printMeshData(obj.grid,cellData);
      obj.outstate.writeVTKfile(obj.vtmBlock,sprintf('Interface_%i',obj.id),...
        obj.grid,time,pointData,cellData,[],[]);
      block = obj.vtmBlock;
    end

    function block = writeVTKfile(obj,fac,time)
      block = obj.writeVTK(fac,time);
    end

    function writeSolution(obj,fac,timeID)
      obj.outstate.results(timeID).time = obj.outstate.timeList(timeID);
      for name = obj.constraintNames
        obj.getConstraint(name).writeSolution(fac,timeID);
      end
    end

    function finalizeOutput(obj)
      obj.outstate.finalize();
      for name = obj.constraintNames
        obj.getConstraint(name).finalizeOutput();
      end
    end

    function stateCurr = getState(obj,varName)
      % get current state variables
      if nargin == 1
        stateCurr = obj.state.get('curr');
      elseif nargin == 2
        stateCurr = obj.state.get('curr',varName);
      else
        error("getState:Number of input arguments must be 1 or 2")
      end
    end

    function stateOld = getStateOld(obj,varName)
      % get last converged state variables
      if nargin == 1
        stateOld = obj.state.get('old');
      elseif nargin == 2
        stateOld = obj.state.get('old',varName);
      else
        error("getState:Number of input arguments must be 1 or 2")
      end
    end


    function stateInit = getStateInit(obj,varName)
      % get initial state variables
      if nargin == 1
        stateInit = obj.state.get('init');
      elseif nargin == 2
        stateInit = obj.state.get('init',varName);
      else
        error("getState:Number of input arguments must be 1 or 2")
      end
    end



    function setState(obj,val,varName)
      % get current state variables
      if nargin == 2
        obj.state.set('curr',val);
      elseif nargin == 3
        obj.state.set('curr',val,varName);
      else
        error("setState:Number of input arguments must be 2 or 3")
      end
    end

    function setStateOld(obj,val,varName)
      % get last converged state variables
      if nargin == 2
        obj.state.set('old',val);
      elseif nargin == 3
        obj.state.set('old',val,varName);
      else
        error("setState:Number of input arguments must be 2 or 3")
      end

    end


    function setStateInit(obj,val,varName)
      % get initial state variables
      if nargin == 2
        obj.state.set('init',val);
      elseif nargin == 3
        obj.state.set('init',val,varName);
      else
        error("setState:Number of input arguments must be 2 or 3")
      end
    end

  end

  methods (Access=private)
    function domains = getInterfaceGrid(obj,in,varargin)

      default = struct('masterSurface',missing,...
        'slaveSurface',missing);

      optionArgs = varargin;
      if isa(in,'Grid')
        optionArgs = varargin(2:end);
      end
      options = readInput(default,optionArgs{:});


      % Process the Interface geometry according to three input formats
      switch class(in)
        case 'FractureGrid'
          % input is a fracture, require the second input to be its hosting
          % domain
          assert(numel(varargin) >= 1 && isa(varargin{1},'Domain') && ...
            isscalar(varargin{1}), 'A fracture requires its hosting Domain');

          matrix = in.matrix;
          neigh = in.cells.matrixFaces;
          faceId = matrix.cells.surfaces.faceId;
          grid1 = matrix.getSurfaceGrid(ismember(faceId,neigh(:,1)));
          grid2 = matrix.getSurfaceGrid(ismember(faceId,neigh(:,2)));
          gridInterf = in;

          domains = [varargin{1};varargin{1}];

        case {'Grid','StructuredGrid'}

          domains = varargin{1};
          assert(isa(domains,'Domain') && numel(domains) == 2, ...
            'An interface grid requires two connected domains');
          domains = domains(:);

          grid1 = getSurfaceGrid(domains(MortarSide.slave).grid,options.slaveSurface);
          grid2 = getSurfaceGrid(domains(MortarSide.master).grid,options.masterSurface);
          gridInterf = in;

        case 'Domain'

          assert(numel(in)==2)

          domains = [in(1);in(2)];

          % the interface grid is made coincide with the slave one!
          grid1 = getSurfaceGrid(domains(MortarSide.slave).grid,options.slaveSurface);
          gridInterf = grid1;
          grid2 = getSurfaceGrid(domains(MortarSide.master).grid,options.masterSurface);

        otherwise
          error('Unexpected input for Interface')

      end


      % process the interface grid supplying all ingredients
      if ~gridInterf.isGeometryProcessed()
        gridInterf.processGeometry();
      end
      obj.grid = InterfaceGrid(gridInterf,grid1,grid2);

      if isempty(obj.grid.neighbors)
        error("No connections found for interface")
      end

    end


  end
end
