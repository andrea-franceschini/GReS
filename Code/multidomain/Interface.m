classdef Interface < handle
  % General Interface class
  %
  % There exist several way in which an interface


  properties (GetAccess=public, SetAccess=public)

    % the grid of the interface (with special information for connection of
    % neighboring grid)
    grid

    % print utilities (pointer to property of solutionScheme)
    outstate

    % interface state
    state

    % id of this interface in the solution scheme object
    id

    % constraints object added to the interface
    constraints

    % name of the constraints 
    constraintNames

    % degree of freedom manager of the interface
    dofm

  end

  properties (Access=private)

    % id of this interface for master and slave domain
    interfaceId

  end


  properties (GetAccess=public, SetAccess=protected)

    % handle to master and slave Discretizer objects
    domains
  end



  methods

    function obj = Interface(input,varargin)

      % domain:  handle to the Discretizer object storing all the
      % information of the model

      % minimal required input is the id of the connected domains

      obj.domains = getInterfaceGrid(obj,input,varargin{:});

      for side = [MortarSide.slave,MortarSide.master]

        obj.domains(side).interfaces = [obj.domains(side).interfaces; obj];

      end

      obj.state = State();


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
        case 'FracturedGrid'
          % input is a fracture, require the second input to be its hosting
          % domain
          assert(isa(varargin{1}),'Domain');

          neigh = in.matrix.cells.matrixFaces;
          faceId = matrix.cells.surfaces.faceId;
          surf = zeros(size(neigh));
          surf(:,1) = find(ismember(faceId,neigh(:,1)));
          surf(:,2) = find(ismember(faceId,neigh(:,2)));
          grid1 = matrix.getSurfaceGrid(surf(:,1));
          grid2 = matrix.getSurfaceGrid(surf(:,2));
          gridInterf = in;

          domains = [varargin{1};varargin{1}];

        case {'Grid','StructuredGrid'}

          domains = [varargin{1}];
          assert(isa(domains,'Domain'))

          grid1 = getSurfaceGrid(domains(MortarSide.slave).grid,options.slaveSurface);
          grid2 = getSurfaceGrid(domains(MortarSide.master).grid,options.masterSurface);
          gridInterf = in;

        case 'Domain'

          assert(numel(in)==2)

          domains = [in(1);in(2)];

          % the interface grid is made coincide with the slave one!
          gridInterf = getSurfaceGrid(domains(MortarSide.slave).grid,options.slaveSurface);
          grid1 = getSurfaceGrid(domains(MortarSide.slave).grid,options.slaveSurface);
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


    function advanceState(obj)

      % note: state in interface solver is just a value struct
      setStateOld(obj,getState(obj));

    end

    function timeStepSetup(obj)
      % prepare the interface solver at the begin of each time step
    end

    function goBackState(obj)

      setState(obj,getStateOld(obj));

    end

    function hasConfigurationChanged = updateConfiguration(obj)

      % base class implements no configuration change
      hasConfigurationChanged = false;

    end

    function initialize(obj)
      % initialize the interface solver

      obj.Jconstraint = [];
      obj.rhsConstraint = [];

      % remove slave Dirichlet boundary conditions for nodal multipliers
      removeSlaveBCents(obj);
    end


    function vtmBlock = writeVTKfile(obj,fac,time)

      grid = obj.grids(MortarSide.slave);
      surfData2D = struct('name', [], 'data', []);
      pointData2D = struct('name', [], 'data', []);
      [surfData,pointData] = writeVTK(obj,fac,time);
      surfData2D = OutState.mergeOutFields(surfData2D,surfData);
      pointData2D = OutState.mergeOutFields(pointData2D,pointData);
      vtmBlock = obj.outstate.vtkFile.createElement('Block');
      obj.outstate.writeVTKfile(vtmBlock,obj.getOutName(),grid,...
        time, pointData2D, surfData2D, [], []);

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




    function J = getJacobian(obj)

      % return the multiplier jacobian block
      J = obj.Jconstraint;

    end

    function rhs = getRhs(obj)

      % return the multiplier jacobian block
      rhs = obj.rhsConstraint;

    end

    function dofs = getMultiplierDoF(obj,el)

      meshSlave = obj.grids(MortarSide.slave);
      dofSlave = getDoFManager(obj,MortarSide.slave);
      ncomp = dofSlave.getNumberOfComponents(obj.coupledVariables);

      if nargin > 1

        dofs = getIncidenceID(obj.multiplierLocation,meshSlave,entityField.cell,el);
        dofs = dofs.getData;

        % dofs = getEntityFromElement( obj.multiplierLocation,...
        %   entityField.cell,...
        %   meshSlave,...
        %   el,...
        %   ncomp );
      else
        dofs = 1:getNumberOfEntities(obj.multiplierLocation,meshSlave);
      end

      dofs = DoFManager.dofExpand(dofs,ncomp);
    end

    function nDoFs = getNumbDoF(obj)

      if ~isempty(obj.nMult)
        nDoFs = obj.nMult;
      else
        ncomp = max(getDoFManager(obj,MortarSide.slave).getNumberOfComponents(obj.coupledVariables));
        obj.nMult = ncomp * getNumberOfEntities(obj.multiplierLocation,...
          obj.grids(MortarSide.slave));
        nDoFs = obj.nMult;
      end
    end


    function R = getRotationMatrix(obj,side,el)

      if obj.multiplierLocation == entityField.cell
        R = obj.grids(side).cells.rotationMatrices(el,:);
        R = reshape(R,3,3);
      elseif obj.multiplierLocation == entityField.node
        nodes = obj.grids(side).getCellNodes(el);
        n = (obj.grids(side).cells.avgNodNormal(nodes,:))';
        n = n(:);
        sz = size(n,1);
        R = zeros(sz(1));
        for i = 1:3:sz(1)
          R(i:i+2,i:i+2) = mxComputeRotationMat(n(i:i+2));
        end
      end
    end


    %%% Get the coupling jacobian
    function J = getJum(obj,side,varargin)

      [s,varId] = getSideAndVar(obj,side,varargin{:});

      J = obj.domains(s).Jum{obj.interfaceId(s)}{varId};
    end

    function J = getJmu(obj,side,varargin)

      [s,varId] = getSideAndVar(obj,side,varargin{:});

      J = obj.domains(s).Jmu{obj.interfaceId(s)}{varId};
    end


    %%% Set the coupling jacobian
    function setJmu(obj,side,setVal,varargin)

      [s,varId] = getSideAndVar(obj,side,varargin{:});

      obj.domains(s).Jmu{obj.interfaceId(s)}{varId} = setVal;
    end

    function setJum(obj,side,setVal,varargin)

      [s,varId] = getSideAndVar(obj,side,varargin{:});

      obj.domains(s).Jum{obj.interfaceId(s)}{varId} = setVal;
    end


    %%% Add contribution to coupling jacobian
    function addJum(obj,side,setVal,varargin)

      [s,varId] = getSideAndVar(obj,side,varargin{:});

      % assert(~isempty(obj.domains(s).Jmu{varId}),...
      %   "Cannot add contribution to empty Jacobian block");

      if ~isempty(obj.domains(s).Jum{obj.interfaceId(s)}{varId})
        obj.domains(s).Jum{obj.interfaceId(s)}{varId} = ...
          obj.domains(s).Jum{obj.interfaceId(s)}{varId} + setVal;
      else
        obj.domains(s).Jum{obj.interfaceId(s)}{varId} = setVal;
      end
    end

    function addJmu(obj,side,setVal,varargin)

      [s,varId] = getSideAndVar(obj,side,varargin{:});

      if ~isempty(obj.domains(s).Jmu{obj.interfaceId(s)}{varId})
        obj.domains(s).Jmu{obj.interfaceId(s)}{varId} = ...
          obj.domains(s).Jmu{obj.interfaceId(s)}{varId} + setVal;
      else
        obj.domains(s).Jmu{obj.interfaceId(s)}{varId} = setVal;
      end


    end


    %%% Add contribution to rhs
    function addRhs(obj,side,setVal,varargin)

      [s,varId] = getSideAndVar(obj,side,varargin{:});

      assert(~isempty(obj.domains(s).rhs{varId}),...
        "Cannot add contribution to empty rhs block");

      obj.domains(s).rhs{varId} = obj.domains(s).rhs{varId} + setVal;
    end

    function var = getVariableField(obj,side,varargin)

      s = getSide(side);

      if nargin < 3
        assert(isscalar(obj.coupledVariables),...
          "The object couples multiple variable field. " + ...
          "Specify the varIable name.")
        var = obj.coupledVariables;
      else
        var = varargin{1};
      end

      dof = getDoFManager(obj,side);
      ents = dof.getActiveEntities(var,1);
      var = getState(obj.domains(s),var);
      var = var(ents);

    end

    function dofm = getDoFManager(obj,side)
      dofm = obj.domains(side).dofm;
    end

    function [out] = isLinear(obj)
       out = false;
    end

  end



  methods (Access=protected)

    function setCoupledVariables(obj,input)

      obj.coupledVariables = obj.getCoupledVariables();

      varMaster = getVariableNames(obj.domains(MortarSide.master).dofm);
      varSlave = getVariableNames(obj.domains(MortarSide.slave).dofm);

      if any([isempty(varMaster),isempty(varSlave)])
        error(['Interface solver can be registered only after ' ...
          'physicsSolver have been defined'])
      end

      % default candidate interface variable are those in common between
      % master and slave discretizers
      sharedVars = intersect(varSlave,varMaster);

      if isempty(obj.coupledVariables)
        in = readInput(struct('variable',sharedVars),input);
        obj.coupledVariables = in.variable;
      end


      isInterfaceValid = all(ismember(obj.coupledVariables,sharedVars));

      if ~isInterfaceValid
        error("The interface attempts to couple a variable named '%s' that is not" + ...
          " available in any of the connected domains.",obj.coupledVariables)
      end

    end

    function removeSlaveBCents(obj)

      nodSlave = obj.grids(MortarSide.slave).cells.loc2glob;

      bc = obj.domains(MortarSide.slave).bcs;
      bcList = bc.getBCList();


      for bcId = bcList
        if getType(bc,bcId) == BCtype.dirichlet
          bcNodes = intersect(bc.getTargetEntities(bcId),nodSlave);
          obj.dirNodes = [obj.dirNodes; bcNodes];
          if obj.multiplierLocation == entityField.node
            bc.removeTargetEntities(bcId,nodSlave);
          end
        end
      end
    end

    function initializeMortarInterface(obj,slaveTag,masterTag)

      surf = [slaveTag,masterTag];

      for side = [MortarSide.slave, MortarSide.master]

        gDom = obj.domains(side).grid;

        if ismissing(surf(side))
          obj.grids(side) = gDom.getSurfaceGrid();
        else
          obj.grids(side) = gDom.getSurfaceGrid(surf(side));
        end
      end

    end


    function checkInterfaceDisjoint(obj)
      % check that the nodes of mortar and slave side are disjoint

      if obj.domainId(1) == obj.domainId(2)
        % interface defined within the same domain

        out = setdiff(obj.interfMesh.local2glob{1},obj.interfMesh.local2glob{2});

        if ~all(ismember(obj.interfMesh.local2glob{1},out))
          error('Nodes of master and slave side are not disjoint');
        end

      end
    end


    function [Ju,Jm] = getCouplingBlocks(obj,id)

      % id: logical array to determine if the request for coupling blocks
      % come from master or slave

      if all(id)
        domId = 1;
      else
        domId = getSide(obj,id);
      end


      if all(id)

        Ju = cell2matrix(obj.Jum{1}) + cell2matrix(obj.Jum{2});
        Jm{1,iVar} = cell2matrix(obj.Jmu{1}) + cell2matrix(obj.Jmu{2});

      else

        Ju{iVar,1} = cell2matrix(obj.Jum{domId});
        Jm{1,iVar} = cell2matrix(obj.Jmu{domId});

      end


    end


    function [s,v] = getSideAndVar(obj,side,varargin)

      s = getSide(side);
      if nargin < 3
        assert(isscalar(obj.coupledVariables),...
          "The object couples multiple variable field. " + ...
          "Specify the varIable name.")
        v = getDoFManager(obj,side).getVariableId(obj.coupledVariables);
      else
        v = getDoFManager(obj,side).getVariableId(varargin{1});
      end

    end

    function outName = getOutName(obj)

      outName = sprintf('Interface_%i',obj.interfId);

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
