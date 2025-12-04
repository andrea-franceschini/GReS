classdef (Abstract) InterfaceSolver < handle
  % Interface to create an interface between domains in GReS

  % Any interfaceSolver should implemented everything needed to couple a
  % variable field across two non conforming lower dimensional (2D) interfaces

  % The interface can belong to two different domains 

  % The fields property lists a set of single physical moduls that are
  % implemented

  % If the physicsSolver is coupled:
  % fields is a string array with all the single physics fields
  % J and rhs are cell variable with the same size as the number of fields


  properties (GetAccess=public, SetAccess=public)
       
    % constraint multiplier block
    Jconstraint

    % rhs block of the constraint equation
    rhsConstraint

    % object for non-conforming integration
    quadrature

    % object holding the interface mesh geometrical info
    interfMesh

    % Number of multiplier dofs
    nMult

    % Location of the lagrange multiplier ( default is P0 )
    multiplierLocation = entityField.surface

    % print utilities
    outstate

    state
    stateOld

    % id of connected domains
    domainId

    % list of variables coupled from the solver
    coupledVariables

    % list of dirichlet nodes - TO DO: consider more elegant handling of dirichlet bcs 
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

    % intialize the state object
    registerInterface(obj);

    % assemble the constraint matrices
    assembleConstraint(obj);

    % update the state after solving a linear system
    updateState(obj,du);

    % update the output structures for printing purposes
    [cellData,pointData] = writeVTK(obj,t);

    % write history to MAT-file
    writeMatFile(obj,t,tID);

  end

  methods (Abstract, Static)

    % return the variables that are coupled from the interface solver
    % if var = empty, the user can/must choose which variable can be coupled
    var = getCoupledVariables(obj)

  end

  methods

    function obj = InterfaceSolver(id,domains,inputStruct)

      % domain:  handle to the Discretizer object storing all the
      % information of the model

      if ~isstruct(inputStruct)
        inputStruct = readstruct(inputStruct,AttributeSuffix="");
      end

      if isfield(inputStruct,"Interface")
        inputStruct = inputStruct.Interface;
      end

      obj.domainId = [getXMLData(inputStruct.Master,[],"domainId");
        getXMLData(inputStruct.Slave,[],"domainId")];

      % store handle to connected domains
      obj.domains = [domains(obj.domainId(1));
        domains(obj.domainId(2))];

      % add handle of this object to the list of interfaces of connected
      % domains
      master = getSide(MortarSide.master);
      slave = getSide(MortarSide.slave);

      obj.interfaceId = [numel(obj.domains(master).interfaces),...
        numel(obj.domains(slave).interfaces)] + 1;

      obj.domains(1).interfaces{obj.interfaceId(master)} = obj;
      obj.domains(2).interfaces{obj.interfaceId(slave)} = obj;
      obj.domains(1).interfaceList(obj.interfaceId(master)) = id;
      obj.domains(2).interfaceList(obj.interfaceId(slave)) = id;


      % initialize interfaces jacobian and rhs blocks to match same number
      % of variable of connected domains
      dofMaster = getDoFManager(obj,MortarSide.master);
      dofSlave = getDoFManager(obj,MortarSide.slave);
      nVarMaster = getNumberOfVariables(dofMaster);
      nVarSlave =  getNumberOfVariables(dofSlave);

      % initialize the jacobian blocks for domain coupling
      obj.domains(master).Jum{obj.interfaceId(master)} = cell(nVarMaster,1);
      obj.domains(slave).Jum{obj.interfaceId(slave)} = cell(nVarSlave,1);

      obj.domains(master).Jmu{obj.interfaceId(master)} = cell(1,nVarMaster);
      obj.domains(slave).Jmu{obj.interfaceId(slave)} = cell(1,nVarSlave);

      % specify the variables to be coupled
      setCoupledVariables(obj,inputStruct)
  
      setMortarInterface(obj,inputStruct);

      % set print utilities to match the slave output (privisonal)
      obj.outstate = OutState(getMesh(obj,MortarSide.slave),inputStruct);

    end


    function advanceState(obj)

      % note: state in interface solver is just a value struct
      % we don't need a copy method to create a deep copy
      obj.stateOld = obj.state;
      
    end

    function goBackState(obj)

      obj.state = obj.stateOld;

    end


    function printState(obj)
      % print solution at the interface according to the print time in the
      % input list 

      if obj.outstate.timeID <= length(obj.outstate.timeList)

        time = obj.outstate.timeList(obj.outstate.timeID);

        % loop over print times within last time step
        while time <= obj.state.t

          assert(time >= obj.stateOld.t, 'Print time %f out of range (%f - %f)',...
            time, obj.stateOld.t, obj.state.t);

          % assert(obj.state.t - obj.stateOld.t > eps('double'),...
          %   'Time step is too small for printing purposes');

          % compute factor to interpolate current and old state variables
          fac = (time - obj.stateOld.t)/(obj.state.t - obj.stateOld.t);
          if isnan(fac) || isinf(fac)
            fac = 1;
          end

          % call methods inside the individual interface solvers
          if obj.outstate.writeVtk
            surfData2D = struct('name', [], 'data', []);
            pointData2D = struct('name', [], 'data', []);
            [surfData,pointData] = writeVTK(obj,fac,time);
            surfData2D = OutState.mergeOutFields(surfData2D,surfData);
            pointData2D = OutState.mergeOutFields(pointData2D,pointData);
            obj.outstate.VTK.writeVTKFile(time, [], [], surfData2D, pointData2D);
          end

          if obj.outstate.writeSolution
            writeMatFile(obj,fac,obj.outstate.timeID);
          end

          obj.outstate.timeID = obj.outstate.timeID + 1;

          if obj.outstate.timeID > length(obj.outstate.timeList)
            break
          else
            time = obj.outstate.timeList(obj.outstate.timeID);
          end

        end

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

      meshSlave = getMesh(obj,MortarSide.slave);
      dofSlave = getDoFManager(obj,MortarSide.slave);
      ncomp = dofSlave.getNumberOfComponents(obj.coupledVariables);

      if nargin > 1

      dofs = getEntityFromElement( obj.multiplierLocation,...
                                   entityField.surface,...
                                   meshSlave,...
                                   el,...
                                   ncomp );
      else
        dofs = 1:getNumberOfEntities(obj.multiplierLocation,meshSlave);
        dofs = DoFManager.dofExpand(dofs,ncomp);
      end
    end

    function nDoFs = getNumbDoF(obj)

      if ~isempty(obj.nMult)
        nDoFs = obj.nMult;
      else
        ncomp = max(getDoFManager(obj,MortarSide.slave).getNumberOfComponents(obj.coupledVariables));
        obj.nMult = ncomp * getNumberOfEntities(obj.multiplierLocation,...
                                                obj.getMesh(MortarSide.slave));
        nDoFs = obj.nMult;
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


    % Set the coupling jacobian
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
      s = getSide(side);
      dofm = obj.domains(s).dofm;
    end

    function mesh = getMesh(obj,side)
      s = getSide(side);
      mesh = obj.interfMesh.msh(s);
    end

  end



  methods (Access=private)

    function setCoupledVariables(obj,input)

      obj.coupledVariables = obj.getCoupledVariables();

      varMaster = getVariableNames(obj.domains(1).dofm);
      varSlave = getVariableNames(obj.domains(2).dofm);
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

    function setMortarInterface(obj,interfaceInput)

      masterSurf = getXMLData(interfaceInput.Master,0,"surfaceTag");
      slaveSurf = getXMLData(interfaceInput.Slave,0,"surfaceTag");

      obj.interfMesh = InterfaceMesh(obj.domains(1).grid.topology,...
        obj.domains(2).grid.topology,...
        masterSurf,slaveSurf);

      checkInterfaceDisjoint(obj);

      multType = getXMLData(interfaceInput.Quadrature,"P0","multiplierType");

      switch multType
        case {"standard","dual"}
          obj.multiplierLocation = entityField.node;
        case "P0"
          obj.multiplierLocation = entityField.surface;
        otherwise
          error("Multiplier type %s is not available." + ...
            "Available types are: standard,dual,P0");
      end

      obj.quadrature = feval(interfaceInput.Quadrature.type,...
                             obj,...
                             multType,...
                             interfaceInput.Quadrature);

      % remove slave Dirichlet boundary conditions for nodal multipliers
      removeSlaveBCents(obj);

      % Mortar quadrature preprocessing
      computeMortarInterpolation(obj);

      % Finalize interface geometry
      finalizeInterfaceGeometry(obj.interfMesh,obj);

      % initialize the state objects
      obj.state.t = 0;
      obj.stateOld = obj.state;
    end


    function computeMortarInterpolation(obj)

      processMortarPairs(obj.quadrature);

      inactiveMaster = ~ismember(1:getMesh(obj,MortarSide.master).nSurfaces,...
        obj.quadrature.interfacePairs(:,2));

      [~, ~, obj.quadrature.interfacePairs(:,2)] = ...
        unique(obj.quadrature.interfacePairs(:,2));

      inactiveSlave = ~ismember(1:getMesh(obj,MortarSide.slave).nSurfaces,...
        obj.quadrature.interfacePairs(:,1));

      [~, ~, obj.quadrature.interfacePairs(:,1)] = ...
        unique(obj.quadrature.interfacePairs(:,1));

      % remove master elements
      obj.interfMesh.removeMortarSurface(1,inactiveMaster);

      % remove slave elements
      obj.interfMesh.removeMortarSurface(2,inactiveSlave);


    end


    function removeSlaveBCents(obj)

      if obj.multiplierLocation ~= entityField.node
        return
      end

      % domain 2 is the slave
      nodSlave = obj.interfMesh.local2glob{2};

      bc = obj.domains(2).bcs;
      bcList = bc.db.keys;

      bcNodes = [];

      for bcId = string(bcList)
        if strcmp(getType(bc,bcId),"Dirichlet")
          bcNodes = [bcNodes; bc.removeBCentities(bcId,nodSlave)];
        end
      end

      obj.dirNodes = bcNodes;

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

      % shortcut to add a contribution to the coupling rhs

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

  end


end
