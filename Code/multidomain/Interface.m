classdef Interface < handle
  % Geometry, state, and system blocks shared by interface constraints.

  properties
    grid
    outstate
    state
    id
    dofm
    J
    rhs
    vtmBlock
  end

  properties (SetAccess=private)
    domains
    constraints
    constraintNames = strings(1,0)
  end

  properties (Access=private)
    interfaceId
    initialized = false
  end

  methods
    function obj = Interface(input,varargin)
      obj.domains = getInterfaceGrid(obj,input,varargin{:});
      obj.dofm = DoFManager(obj.grid);
      obj.state = State();
      obj.constraints = containers.Map('KeyType','char','ValueType','any');

      for side = [MortarSide.slave,MortarSide.master]
        domain = obj.domains(side);
        slot = numel(domain.interfaces) + 1;
        domain.interfaces{slot} = obj;
        obj.interfaceId(side) = slot;
      end
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

    function initialize(obj)
      nFields = obj.dofm.getNumberOfVariables();
      obj.J = cell(nFields);
      obj.rhs = cell(nFields,1);

      for side = [MortarSide.slave,MortarSide.master]
        [domain,slot] = obj.getDomainSlot(side);
        nDomainFields = domain.dofm.getNumberOfVariables();
        domain.Jum{slot} = cell(nDomainFields,nFields);
        domain.Jmu{slot} = cell(nFields,nDomainFields);
      end

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

    function J = getJacobian(obj,varargin)
      narginchk(1,3);
      J = obj.J;
      if ~isempty(varargin)
        rows = obj.getFieldIds(obj.dofm,varargin{1});
        cols = rows;
        if numel(varargin) == 2
          cols = obj.getFieldIds(obj.dofm,varargin{2});
        end
        J = J(rows,cols);
      end
    end

    function rhs = getRhs(obj,varargin)
      narginchk(1,2);
      rhs = obj.rhs;
      if ~isempty(varargin)
        rhs = rhs(obj.getFieldIds(obj.dofm,varargin{1}));
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

    function J = getJum(obj,side,domainField,interfaceField)
      [domain,slot] = obj.getDomainSlot(side);
      J = domain.Jum{slot};
      if nargin >= 3
        J = J(obj.getFieldIds(domain.dofm,domainField),:);
      end
      if nargin == 4
        J = J(:,obj.getFieldIds(obj.dofm,interfaceField));
      end
      if isscalar(J) && nargin >= 3
        J = J{1};
      end
    end

    function J = getJmu(obj,side,domainField,interfaceField)
      [domain,slot] = obj.getDomainSlot(side);
      J = domain.Jmu{slot};
      if nargin >= 3
        J = J(:,obj.getFieldIds(domain.dofm,domainField));
      end
      if nargin == 4
        J = J(obj.getFieldIds(obj.dofm,interfaceField),:);
      end
      if isscalar(J) && nargin >= 3
        J = J{1};
      end
    end

    function setJum(obj,side,val,varargin)
      [domain,slot,r,c] = obj.getCouplingIndices(side,varargin{:});
      domain.Jum{slot}{r,c} = val;
    end

    function setJmu(obj,side,val,varargin)
      [domain,slot,r,c] = obj.getCouplingIndices(side,varargin{:});
      domain.Jmu{slot}{c,r} = val;
    end

    function addJum(obj,side,val,varargin)
      [domain,slot,r,c] = obj.getCouplingIndices(side,varargin{:});
      domain.Jum{slot}{r,c} = obj.addBlock(domain.Jum{slot}{r,c},val);
    end

    function addJmu(obj,side,val,varargin)
      [domain,slot,r,c] = obj.getCouplingIndices(side,varargin{:});
      domain.Jmu{slot}{c,r} = obj.addBlock(domain.Jmu{slot}{c,r},val);
    end

    function addRhs(obj,side,val,domainField)
      domain = obj.domains(side);
      id = domain.dofm.getVariableId(domainField);
      domain.rhs{id} = obj.addBlock(domain.rhs{id},val);
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
    function [domain,slot] = getDomainSlot(obj,side)
      domain = obj.domains(side);
      slot = obj.interfaceId(side);
    end

    function [domain,slot,r,c] = getCouplingIndices(obj,side,domainField,interfaceField)
      [domain,slot] = obj.getDomainSlot(side);
      if nargin < 4
        assert(obj.dofm.getNumberOfVariables() == 1, ...
          "Specify the interface field when there are multiple interface fields");
        interfaceField = 1;
      end
      r = domain.dofm.getVariableId(domainField);
      c = obj.dofm.getVariableId(interfaceField);
    end

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

  methods (Static, Access=private)
    function ids = getFieldIds(dofm,fields)
      if isnumeric(fields)
        ids = dofm.getVariableId(fields);
      else
        ids = arrayfun(@(field) dofm.getVariableId(field),string(fields));
      end
    end

    function block = addBlock(block,val)
      if isempty(block)
        block = val;
      else
        block = block + val;
      end
    end
  end

  methods (Static)

    function interfaces = addInterfaces(domains,input)

      assert(nargin == 2,"Input must be a structure array or an input file")

      interfStruct = readInput(input);

      interfaces = {};

      interfNames = fieldnames(interfStruct);

      for i = 1:numel(interfNames)

        % deal with multiple interfaces having same name
        for in = [interfStruct.(interfNames{i})]

          interfaces = Interface.add(interfNames{i},...
            domains,...
            interfaces,in);
        end
      end

    end


    function interfaces = add(interfType,domains,varargin)
      %   Add a new interface object to the interface list.
      %
      %   interfaces = add(interfType, domains, input)
      %   interfaces = add(interfType, domains, interfaces, input)
      %
      %   This function creates a new interface object of type `interfType`,
      %   assigns it a progressive index, registers solver-specific data,
      %   and appends it to the interface list.
      %
      %   INPUT
      %   -----
      %   interfType : function handle or class name (char/string)
      %       Constructor for the interface object. The constructor must
      %       accept the signature:
      %           obj = interfType(id, domains, input)
      %
      %   domains : array or cell array
      %       Collection of domains connected by the interface.
      %
      %   interfaces : cell array (optional)
      %       Existing list of interface objects. If omitted, a new list
      %       is created.
      %
      %   input : cell array
      %       Additional user-defined input parameters passed both to the
      %       constructor and to the interface method `registerInterface`.
      %
      %   OUTPUT
      %   ------
      %   interfaces : cell array
      %       Updated list of interface objects, including the newly
      %       created interface appended at the end.
      %
      %
      %   See also: Interface
      %
      %   ---------------------------------------------------------------------


      if nargin == 2
        interfaces = {};
        input = [];
      elseif iscell(varargin{1})
        interfaces = varargin{1};
        input = varargin{2:end};
      else
        % generate the interface list
        interfaces = {};
        input = varargin{1:end};
      end

      k = numel(interfaces);



      % general input
      interf = feval(interfType,k+1,domains,input);
      % solver-specific input
      interf.registerInterface(input);

      interfaces{end+1} = interf;

    end
  end
end
