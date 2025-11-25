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

  end

  properties (Access=private)

    % id of connected domains
    domainId

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

    % advance the state after reaching time step convergence
    advanceState(obj);

  end

  methods

    function obj = InterfaceSolver(domains,inputStruct)

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
      obj.interfaceId = [numel(obj.domains(1).interfaces),...
                         numel(obj.domains(2).interfaces)] + 1;

      obj.domains(1).interfaces{obj.interfaceId(1)} = obj;
      obj.domains(2).interfaces{obj.interfaceId(2)} = obj;

      % initialize interfaces jacobian and rhs blocks to match same number
      % of variable of connected domains
      nVarMaster = getNumberOfVariables(obj.domains(1).dofm);
      nVarSlave =  getNumberOfVariables(obj.domains(2).dofm);

      % initialize the jacobian blocks for domain coupling
      obj.domains(1).Jum{end+1} = cell(nVarMaster,1);
      obj.domains(2).Jum{end+1} = cell(nVarSlave,1);

      obj.domains(1).Jmu{end+1} = cell(1,nVarMaster);
      obj.domains(2).Jmu{end+1} = cell(1,nVarSlave);

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

      multType = getXMLData(interfaceInput,"P0","multiplierType");

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



    function J = getJacobian(obj)

      % return the multiplier jacobian block
      J = obj.Jconstraint;

    end

    function dofs = getMultiplierDoF(obj,el)

      ncomp = obj.domains(2).dofm.getNumberOfComponents(obj.coupledVariable);

      dofs = getEntityFromElement( obj.multiplierLocation,...
                                   entityField.surface,...
                                   obj.interfMesh.msh(2),...
                                   el,...
                                   ncomp);
    end

    function nDoFs = getNumbDoF(obj)

      nDoFs = obj.nMult;
    end


    % set methods to ease the access to the domain coupling matrices

    function setJum(obj,side,setVal,varargin)

      s = getSide(side);
      if nargin < 4
        assert(isscalar(obj.coupledVariable),...
          "The objace couples multiple variable field. " + ...
          "Specify the varIable name.")
        varId = getDoF(side).getVariableId(obj.coupledVariable);
      else
        varId = getDoF(side).getVariableId(varargin{1});
      end

      obj.domains(s).Jum{obj.interfaceId(s)}{varId} = setVal;
    end


    function setJmu(obj,side,setVal,varargin)

      s = getSide(side);
      if nargin < 4
        assert(isscalar(obj.coupledVariable),...
          "The objace couples multiple variable field. " + ...
          "Specify the varIable name.")
        varId = getDoF(side).getVariableId(obj.coupledVariable);
      else
        varId = getDoF(side).getVariableId(varargin{1});
      end

      obj.domains(s).Jmu{obj.interfaceId(s)}{varId} = setVal;

    end

    function addRhs(obj,side,setVal,varargin)

      % shortcut to add a contribution to the coupling rhs 

      s = getSide(side);
      if nargin < 4
        assert(isscalar(obj.coupledVariable),...
          "The objace couples multiple variable field. " + ...
          "Specify the varIable name.")
        varId = getDoF(side).getVariableId(obj.coupledVariable);
      else
        varId = getDoF(side).getVariableId(varargin{1});
      end

      assert(~isempty(obj.domains(s).rhs{varId}),...
        "Cannot add contribution to empty rhs block");

      obj.domains(s).rhs{varId} = obj.domains(s).rhs{varId} + setVal;

    end


    function dofm = getDoF(side,interf)
      s = getSide(side);
      dofm = interf.domains(s).dofm;
    end

  end



  methods (Access=private)

    
    function computeMortarInterpolation(obj)

      processMortarPairs(obj.quadrature); 

      inactiveMaster = ~ismember(1:obj.interfMesh.msh(1).nSurfaces,...
        obj.quadrature.interfacePairs(:,2));

      [~, ~, obj.quadrature.interfacePairs(:,2)] = ...
        unique(obj.quadrature.interfacePairs(:,2));

      inactiveSlave = ~ismember(1:obj.interfMesh.msh(2).nSurfaces,...
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

      bc = obj.solvers(2).bcs;
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


    function side = getSide(obj,domId)

      side = find(obj.domainId == domId);

      if isempty(side)
        error("Domain %i is not connected by the interface")
      end

    end

  end


end
