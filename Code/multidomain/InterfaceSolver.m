classdef (Abstract) InterfaceSolver < handle
  % General Interface solver between domains in GReS
  % 
  % Any interfaceSolver should implemented everything needed to couple a
  % variable field across two non conforming lower dimensional (2D) interfaces
  %
  % Each side of the interface is specified by one (or more) surfaceTags of
  % a domain. The sides of the interface can both belong to the same
  % domain, provided that their nodes are disjoint.


  properties (GetAccess=public, SetAccess=public)
       
    % constraint multiplier block
    Jconstraint

    % rhs block of the constraint equation
    rhsConstraint

    % object for non-conforming integration
    quadrature

    % grids(1): master    grids(2): slave
    grids

    % Number of multiplier dofs
    nMult

    % Location of the lagrange multiplier ( default is P0 )
    multiplierLocation = entityField.surface

    % print utilities (pointer to property of solutionScheme)
    outstate

    % interface state
    state 
    stateOld

    % id of this interface in the solution scheme object
    interfId

    % global id of the connected domains
    domainId

    % list of variables coupled from the solver
    coupledVariables

    % list of dirichlet nodes
    dirNodes

  end

  properties (Access=private)

    % id of this interface for master and slave domain
    interfaceId

  end


  properties (GetAccess=public, SetAccess=protected)

    % handle to master and slave Discretizer objects
    domains
  end

  methods (Abstract)

    % intialize the interface properties
    registerInterface(obj);

    % assemble the constraint matrices
    assembleConstraint(obj,varargin);

    % update the state after solving a linear system
    updateState(obj,du);

    % update the output structures for printing purposes
    [cellData,pointData] = writeVTK(obj,t);

    % write history to MAT-file
    writeSolution(obj,t,tID);

  end

  methods (Abstract, Static)

    % return the variables that are coupled from the interface solver
    % if var = empty, the user can/must choose which variable can be coupled
    var = getCoupledVariables(obj)

  end

  methods

    function obj = InterfaceSolver(id,domains,varargin)

      % domain:  handle to the Discretizer object storing all the
      % information of the model

      % minimal required input is the id of the connected domains
      default = struct('slaveDomain',[], ...
        'masterDomain',[]);

      input = readInput(default,varargin{:});


      obj.domainId = [input.slaveDomain,input.masterDomain];


      % store handle to connected domains
      obj.domains = [domains(obj.domainId(1));
        domains(obj.domainId(2))];

      % add handle of this object to the list of interfaces of connected
      % domains

      obj.interfaceId(MortarSide.slave) = numel(obj.domains(MortarSide.slave).interfaces) + 1;
      obj.interfaceId(MortarSide.master) = numel(obj.domains(MortarSide.master).interfaces) + 1;


      for side = [MortarSide.slave,MortarSide.master]

        obj.domains(side).interfaces{obj.interfaceId(side)} = obj;
        obj.domains(side).interfaceList(obj.interfaceId(side)) = id;

        dof = getDoFManager(obj,side);
        nVar = getNumberOfVariables(dof);
        obj.domains(side).Jum{obj.interfaceId(side)} = cell(nVar,1);
        obj.domains(side).Jmu{obj.interfaceId(side)} = cell(1,nVar);

      end

      % prepare cross-grid informations
      setMortarInterface(obj,input);

      % specify the variables to be coupled
      setCoupledVariables(obj,input)


    end

    function advanceState(obj)

      % note: state in interface solver is just a value struct
      obj.stateOld = obj.state;
      
    end

    function timeStepSetup(obj)
      % prepare the interface solver at the begin of each time step
    end

    function goBackState(obj)

      obj.state = obj.stateOld;

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
      nList = getFlatConnectivity(grid,"surfaces");
      obj.outstate.writeVTKfile(vtmBlock,obj.getOutName(),nList,...
        time, [], [], pointData2D, surfData2D);

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

        dofs = getIncidenceID(obj.multiplierLocation,meshSlave,entityField.surface,el);
        dofs = dofs.getData;

        % dofs = getEntityFromElement( obj.multiplierLocation,...
        %   entityField.surface,...
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
        R = obj.grid(side).sirfaces.rotationMatrices(elemId,:);
        R = reshape(R,3,3);
      elseif obj.multiplierLocation == entityField.node
        nodes = obj.grids(side).getSurfNodes(el);
        n = (obj.grids(side).surfaces.avgNodNormal(nodes,:))';
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

    function slaveArea = getSlaveArea(obj)
      slaveArea = obj.quadrature.areaSlave;
    end


    function dofm = getDoFManager(obj,side)
      dofm = obj.domains(side).dofm;
    end



  end



  methods (Access=protected)

    function setCoupledVariables(obj,input)

      obj.coupledVariables = obj.getCoupledVariables();

      varMaster = getVariableNames(obj.domains(1).dofm);
      varSlave = getVariableNames(obj.domains(2).dofm);

      if any([isempty(varMaster),isempty(varSlave)])
        error(['Interface solver can be registered only after ' ...
          'physicsSolver have been defined'])
      end

      % default candidate interface variable are those in common between
      % master and slave discretizers
      sharedVars = intersect(varSlave,varMaster);

      if ~isempty(obj.coupledVariables)
        % the interfaceSolver specifies the coupled variables directly in
        % the properties block
        % only check that the interface is compatible with the available
        % field

        isInterfaceValid = all(ismember(obj.coupledVariables,sharedVars));

        if ~isInterfaceValid
          error("The interface attempts to couple a variable that is not" + ...
            " available in any of the connected domains.")
        end

      else
        % the interfaceSolver does not specify the coupled variables
        % in the properties block

        if isfield(input,"variable")
          sharedVars = getXMLData(input,[],"variable");
        end

        obj.coupledVariables = sharedVars;
      end
    end

    function setMortarInterface(obj,params)


      default = struct('masterSurface',missing,...
                       'slaveSurface',missing,...
                       'multiplierType',"P0",...
                       'Quadrature',struct());

      % the Quadrature field implies that this is an xml field

      params = readInput(default,params);

      processMortarGrid(obj,params);

      % initialize the state objects
      obj.state.t = 0;
      obj.stateOld = obj.state;
    end


    function processMortarGrid(obj,params)

      obj.grids = repmat(Grid.empty,2,1);

      % read parameters 
      switch params.multiplierType
        case {"standard","dual"}
          obj.multiplierLocation = entityField.node;
        case "P0"
          obj.multiplierLocation = entityField.surface;
        otherwise
          error("Multiplier type %s is not available." + ...
            "Available types are: standard,dual,P0");
      end

      if numel(fieldnames(params.Quadrature)) > 0
        quadType = params.Quadrature.type;
      else
        quadType = "SegmentBasedQuadrature";
      end


      % get initial grids (may include non active mortar elements)
      tags = {params.slaveSurface,params.masterSurface};

      for side = [MortarSide.slave,MortarSide.master]

        grid = obj.domains(side).grid;
        if ismissing(tags{side})
          tags{side} = (1:grid.surface.nTag)';
        end
        obj.grids(side) = getSurfaceGrid(grid,tags{side});
      end

      cs = ContactSearching(obj.grids(MortarSide.slave),obj.grids(MortarSide.master));

      elemConnectivity = cs.getElementConnectivity();


      obj.quadrature = feval(quadType,...
                             params.multiplierType, ...
                             obj.grids, ...
                             params.Quadrature);


      processMortarPairs(obj.quadrature,elemConnectivity);


      % now finalize the grids and compute edges
      for side = [MortarSide.slave,MortarSide.master]
        
        surf = obj.domains(side).grid.surfaces;
        isMortarElem = false(surf.num,1);
        surfList = find(ismember(surf.tag,tags{side}));
        id = surfList(obj.quadrature.interfacePairs(:,side));
        isMortarElem(id) = true;

        % update the grids with only active mortar elements
        obj.grids(side) = getSurfaceGrid(obj.domains(side).grid,isMortarElem);

        % process the edges of the final grid
        obj.grids(side).processEdges();
        obj.grids(side).computeAvgNodalNormal();
        surf = obj.grids(side).surfaces;
        R = zeros(surf.num,9);
        for i = 1:surf.num
          n = surf.normal(i,:);
          Ri = mxComputeRotationMat(n);
          R(i,:) = Ri(:);
        end
        obj.grids(side).surfaces.rotationMatrices = R;
      end
      
    end



    function removeSlaveBCents(obj)

      nodSlave = obj.grids(MortarSide.slave).surfaces.loc2glob;

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

          interfaces = InterfaceSolver.add(interfNames{i},...
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
      %   See also: InterfaceSolver
      %
      %   ---------------------------------------------------------------------


      if iscell(varargin{1})
        interfaces = varargin{1};
        i1 = 2;
      else
        % generate the interface list
        interfaces = {};
        i1 = 1;
      end

      k = numel(interfaces);

      input = varargin{i1:end};

      % general input
      interf = feval(interfType,k+1,domains,input);
      % solver-specific input
      interf.registerInterface(input);

      interfaces{end+1} = interf;

    end










    end


end
