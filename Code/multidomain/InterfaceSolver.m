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

    % id of connected domains
    domainId


  end

  properties (Access=private)

    % id of this interface for master and slave domain
    interfaceId

  end

  properties (Abstract)

    coupledVariable 
    state
    stateOld

  end


  properties (GetAccess=public, SetAccess=protected)

    % handle to master and slave Discretizer objects
    domains
  end

  methods (Abstract)

    % intialize the state object
    initializeConstraint(obj);

    % assemble the constraint matrices
    assembleConstraint(obj);

    % update the state after solving a linear system
    updateState(obj,du);

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
      obj.domains(master).Jum{end+1} = cell(nVarMaster,1);
      obj.domains(slave).Jum{end+1} = cell(nVarSlave,1);

      obj.domains(master).Jmu{end+1} = cell(1,nVarMaster);
      obj.domains(slave).Jmu{end+1} = cell(1,nVarSlave);

    end

    function setInterface(obj,interfaceInput)

      % setup the interfaceMesh object and process the mortar quadrature

      % the following operations are done in every interface solver

      masterSurf = getXMLData(interfaceInput.Master,0,"surfaceTag");
      slaveSurf = getXMLData(interfaceInput.Slave,0,"surfaceTag");

      obj.interfMesh = interfaceMesh(obj.domains(1).grid.topology,...
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

      % initialize the state object and the number of multipliers
      initializeConstraint(obj)

      % Finalize interface geometry
      finalizeInterfaceGeometry(obj.interfMesh,obj);

      % set print utilities to match the slave output (privisonal)
      obj.outstate = OutState(obj.interfMesh.msh(2),interfaceInput);
      
    end


    function advanceState(obj)

      obj.stateOld = obj.state;
      
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
      ncomp = dofSlave.getNumberOfComponents(obj.coupledVariable);

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

      nDoFs = obj.nMult;
    end


    % set methods to ease the access to the domain coupling matrices


    function J = getJum(obj,side,varargin)

      s = getSide(side);

      if nargin < 3
        assert(isscalar(obj.coupledVariable),...
          "The objace couples multiple variable field. " + ...
          "Specify the varIable name.")
        varId = getDoFManager(obj,side).getVariableId(obj.coupledVariable);
      else
        varId = getDoFManager(obj,side).getVariableId(varargin{1});
      end

      J = obj.domains(s).Jum{obj.interfaceId(s)}{varId};
    end

    function J = getJmu(obj,side,varargin)

      s = getSide(side);

      if nargin < 3
        assert(isscalar(obj.coupledVariable),...
          "The objace couples multiple variable field. " + ...
          "Specify the varIable name.")
        varId = getDoFManager(obj,side).getVariableId(obj.coupledVariable);
      else
        varId = getDoFManager(obj,side).getVariableId(varargin{1});
      end

      J = obj.domains(s).Jmu{obj.interfaceId(s)}{varId};
    end


    function setJmu(obj,side,setVal,varargin)

      s = getSide(side);
      if nargin < 4
        assert(isscalar(obj.coupledVariable),...
          "The objace couples multiple variable field. " + ...
          "Specify the varIable name.")
        varId = getDoFManager(obj,side).getVariableId(obj.coupledVariable);
      else
        varId = getDoFManager(obj,side).getVariableId(varargin{1});
      end

      obj.domains(s).Jmu{obj.interfaceId(s)}{varId} = setVal;

    end

    function setJum(obj,side,setVal,varargin)

      s = getSide(side);
      if nargin < 4
        assert(isscalar(obj.coupledVariable),...
          "The object couples multiple variable field. " + ...
          "Specify the varIable name.")
        varId = getDoFManager(obj,side).getVariableId(obj.coupledVariable);
      else
        varId = getDoFManager(obj,side).getVariableId(varargin{1});
      end

      obj.domains(s).Jum{obj.interfaceId(s)}{varId} = setVal;
    end


    function addRhs(obj,side,setVal,varargin)

      % shortcut to add a contribution to the coupling rhs 

      s = getSide(side);
      if nargin < 4
        assert(isscalar(obj.coupledVariable),...
          "The object couples multiple variable field. " + ...
          "Specify the varIable name.")
        varId = getDoFManager(obj,side).getVariableId(obj.coupledVariable);
      else
        varId = getDoFManager(obj,side).getVariableId(varargin{1});
      end

      assert(~isempty(obj.domains(s).rhs{varId}),...
        "Cannot add contribution to empty rhs block");

      obj.domains(s).rhs{varId} = obj.domains(s).rhs{varId} + setVal;

    end

    function var = getVariableField(obj,side,varargin)

      s = getSide(side);

      if nargin < 3
        assert(isscalar(obj.coupledVariable),...
          "The object couples multiple variable field. " + ...
          "Specify the varIable name.")
        var = obj.coupledVariable;
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

    
    function computeMortarInterpolation(obj)

      processMortarPairs(obj.quadrature); 

      inactiveMaster = ~ismember(1:getMesh(obj,MortarSide.master).nSurfaces,...
        obj.quadrature.interfacePairs(:,2));

      [~, ~, obj.quadrature.interfacePairs(:,2)] = ...
        unique(obj.quadrature.interfacePairs(:,2));

      inactiveSlave = ~ismember(1:getMesh(obj,MortarSide.master).nSurfaces,...
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

      for bcId = string(bcList)
        if strcmp(getType(bc,bcId),"Dirichlet")
          bc.removeBCentities(bcId,nodSlave);
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

  end


end
